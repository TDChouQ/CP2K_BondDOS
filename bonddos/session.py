"""High-level k-point bonding session used by both the GUI and the command line."""
from __future__ import annotations

import csv
import datetime as _dt
import hashlib
import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import __version__
from .bonding import (BondingEngine, BondingOptions, BondResult, gauge_check, resolve_wrap_ambiguity,
                      QLABEL, default_grid)
from .constants import FWHM_PER_SIGMA
from .kspace import KSystem, load_system, Bond
from .population import populations, gross_populations, kpoint_pdos, PDOSGroup
from .spectra import make_grid
from . import lobster


@dataclass
class BondRequest:
    a: int
    b: int
    r: tuple | None = None
    all_images: bool = False
    mode: str = "IMAGE"
    resolve: bool = False

    @staticmethod
    def parse(text: str, default_mode="IMAGE"):
        """'1 2', '1 2 ALL', '1 2 -1 0 0', '1 2 PAIR', optional trailing 'ORB' to orbital-resolve."""
        toks = text.replace(",", " ").split()
        if len(toks) < 2:
            raise ValueError(f"键定义至少需要两个原子编号: {text!r}")
        a, b = int(toks[0]), int(toks[1])
        rest = [t.upper() for t in toks[2:]]
        resolve = False
        if "ORB" in rest:
            resolve = True
            rest.remove("ORB")
        mode = default_mode
        r = None
        allim = False
        if rest and rest[0] == "PAIR":
            mode = "PAIR"
            rest = rest[1:]
        elif rest and rest[0] == "ALL":
            allim = True
            rest = rest[1:]
        if rest:
            if len(rest) != 3:
                raise ValueError(f"晶格平移需要 3 个整数: {text!r}")
            r = tuple(int(x) for x in rest)
            mode = "IMAGE"
        return BondRequest(a, b, r, allim, mode, resolve)


@dataclass
class SessionSettings:
    basis: str = "ao"
    quantities: tuple = ("cohp", "coop", "cobi")
    resolution: str = "l"
    spins: str = "auto"
    fwhm: float = 0.2355            # eV (sigma 0.1)
    shape: str = "gaussian"
    emin: float = -15.0
    emax: float = 10.0
    de: float = 0.01
    eref_mode: str = "auto"         # auto|hoco|midgap|fermi|manual
    eref_value: float | None = None


class BondSession:
    def __init__(self):
        self.system: KSystem | None = None
        self.settings = SessionSettings()
        self.results: list[BondResult] = []
        self.engine: BondingEngine | None = None
        self.eref = 0.0
        self.eref_note = ""
        self.gauge = {}
        self.notes = []
        self.inputs = {}
        self.pop_m = None
        self.pop_l = None

    # ---------------------------------------------------------- loading
    def load(self, mokp, ks_dir=None, ks_set=None, triangle="FULL", periodic="XYZ", electrons=None, progress=None):
        self.system = load_system(mokp, ks_dir or None, ks_set, triangle, periodic,
                                  progress=(lambda f, t: progress(0.8 * f, t)) if progress else None)
        self.inputs = {"mokp": str(mokp), "ks_dir": str(ks_dir) if ks_dir else None, "ks_set": ks_set,
                       "triangle": triangle, "periodic": periodic}
        self.system.validate(progress=(lambda f, t: progress(0.8 + 0.15 * f, t)) if progress else None,
                             electrons=electrons)
        self.notes = resolve_wrap_ambiguity(self.system)
        self.gauge = gauge_check(self.system)
        if progress:
            progress(1.0, "数据加载完成")
        self.update_reference()
        self.results = []
        self.pop_m = self.pop_l = None
        return self.system

    def update_reference(self):
        s = self.settings
        self.eref, self.eref_note = self.system.reference_energy(s.eref_mode, s.eref_value)
        return self.eref

    def summary(self):
        sysm = self.system
        d = sysm.d
        vbm, cbm, partial = sysm.band_edges()
        elements = {}
        for a in d.atoms:
            elements[a.element] = elements.get(a.element, 0) + 1
        return {
            "atoms": len(d.atoms), "formula": " ".join(f"{k}{v}" for k, v in elements.items()),
            "nao": d.nao, "nmo": d.nmo, "nkp": d.nkp, "nspins": d.nspins,
            "mesh": "×".join(str(x) for x in sysm.mesh.mesh), "full_mesh": sysm.mesh.full,
            "mesh_messages": sysm.mesh.messages, "H": sysm.h_source, "S": sysm.has_S,
            "vbm": vbm, "cbm": cbm, "gap": (cbm - vbm) if (vbm is not None and cbm is not None) else None,
            "partial": partial, "electrons": d.electrons(), "eref": self.eref, "eref_note": self.eref_note,
            "use_real_wfn": d.use_real_wfn, "warnings": sysm.warnings + self.notes,
        }

    # ---------------------------------------------------------- bonds
    def build_bonds(self, requests: list[BondRequest], auto=None):
        """auto: dict(dmin, dmax, elements_a, elements_b, atoms) for neighbour search."""
        st = self.system.structure
        eng = BondingEngine(self.system, BondingOptions(self.settings.basis, self.settings.quantities,
                                                         self.settings.resolution, self.settings.spins))
        bonds, resolve = [], set()
        seen = set()
        for rq in requests:
            for b in eng.make_bond(rq.a, rq.b, rq.r, rq.mode, rq.all_images):
                if b.key() in seen:
                    continue
                seen.add(b.key())
                if rq.resolve:
                    resolve.add(len(bonds))
                bonds.append(b)
        if auto:
            for a, b, R, d in st.neighbor_search(auto.get("dmin", 0.0), auto["dmax"], auto.get("elements_a"),
                                                  auto.get("elements_b"), auto.get("atoms")):
                bb = Bond(a, b, tuple(R), st.wrap_translation(a, b, R), d, "IMAGE",
                          f"{st.label(a)}–{st.label(b)} [{R[0]} {R[1]} {R[2]}]")
                if bb.key() in seen:
                    continue
                seen.add(bb.key())
                if auto.get("resolve"):
                    resolve.add(len(bonds))
                bonds.append(bb)
        return eng, bonds, resolve

    def compute(self, requests, auto=None, progress=None, cancel=None):
        self.update_reference()
        eng, bonds, resolve = self.build_bonds(requests, auto)
        if not bonds:
            raise ValueError("没有需要计算的键：请添加键或设置自动搜索的距离范围")
        if len(bonds) > 3000:
            raise ValueError(f"键数量过多 ({len(bonds)})，请缩小距离范围或限定元素")
        self.engine = eng
        self.results = eng.run(bonds, resolve, progress=progress, cancel=cancel)
        return self.results

    # ---------------------------------------------------------- curves
    def grid(self):
        s = self.settings
        return make_grid(s.emin, s.emax, s.de)

    def curves(self, result: BondResult, q, channel=None, cumulative=False, grid=None):
        s = self.settings
        g = self.grid() if grid is None else grid
        return g, self.engine.spectra(result, q, self.eref, g, s.fwhm, s.shape, channel, cumulative)

    def sum_curves(self, results, q, cumulative=False):
        g = self.grid()
        acc = {}
        for r in results:
            _, c = self.curves(r, q, cumulative=cumulative, grid=g)
            for sp, y in c.items():
                acc[sp] = acc.get(sp, 0) + y
        return g, acc

    # ---------------------------------------------------------- populations / pdos
    def compute_populations(self, progress=None):
        if not self.system.has_S:
            raise ValueError("需要 S(k)（MO_KP AO_EXPORT_TYPE OVERLAP_MATRIX）")
        self.pop_m = populations(self.system, "mulliken", progress)
        self.pop_l = populations(self.system, "lowdin", progress)
        return self.pop_m, self.pop_l

    def pdos(self, groups, method="mulliken", include_tdos=True, occupied=False, progress=None):
        s = self.settings
        g = self.grid()
        return g, kpoint_pdos(self.system, groups, self.eref, g, s.fwhm, s.shape, method, include_tdos,
                              progress, occupied)

    # ---------------------------------------------------------- export
    def export(self, folder, prefix="CP2K", lobster_files=True, pdos=True, progress=None):
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        sysm, st = self.system, self.system.structure
        written = []
        spins = self.engine.spins if self.engine else list(range(sysm.ns))
        spin_names = {0: "up", 1: "down"} if sysm.ns == 2 else {0: "total"}
        g = self.grid()
        qs = self.engine.quantities if self.engine else []
        # bond summary
        if self.results:
            p = folder / f"{prefix}_bonds.csv"
            with open(p, "w", newline="", encoding="utf-8-sig") as fh:
                w = csv.writer(fh)
                head = ["No", "label", "atomA", "atomB", "elemA", "elemB", "R1", "R2", "R3", "length_A", "mode", "status"]
                for q in qs:
                    for sp in spins:
                        head.append(f"I{QLABEL[q]}_{spin_names[sp]}")
                    if sysm.ns == 2:
                        head.append(f"I{QLABEL[q]}_total")
                    if q == "cohp":
                        head.append("-ICOHP_total_eV")
                head.append("notes")
                w.writerow(head)
                for i, r in enumerate(self.results, 1):
                    b = r.bond
                    row = [i, b.label, b.a, b.b, st.symbol(b.a), st.symbol(b.b), *b.r_raw, f"{b.length:.6f}", b.mode, r.status]
                    for q in qs:
                        I = r.integrals.get(q, {})
                        for sp in spins:
                            key = spin_names[sp]
                            row.append(f"{I.get(key, float('nan')):.6f}")
                        if sysm.ns == 2:
                            row.append(f"{I.get('total', float('nan')):.6f}")
                        if q == "cohp":
                            row.append(f"{-I.get('total', float('nan')):.6f}")
                    row.append("; ".join(r.notes))
                    w.writerow(row)
            written.append(p)
            # curves
            for q in qs:
                p = folder / f"{prefix}_{QLABEL[q]}_curves.csv"
                cols = [g]
                head = ["E-Eref_eV"]
                for i, r in enumerate(self.results, 1):
                    if r.status in ("aliased", "nonperiodic"):
                        continue
                    _, c = self.curves(r, q, grid=g)
                    _, ci = self.curves(r, q, cumulative=True, grid=g)
                    for sp in spins:
                        tag = f"{i}:{r.bond.label}:{spin_names[sp]}"
                        cols += [c[sp], ci[sp]]
                        head += [f"{QLABEL[q]}[{tag}]", f"I{QLABEL[q]}[{tag}]"]
                with open(p, "w", newline="", encoding="utf-8-sig") as fh:
                    w = csv.writer(fh)
                    w.writerow(head)
                    for row in np.column_stack(cols):
                        w.writerow([f"{x:.6g}" for x in row])
                written.append(p)
            # orbital-resolved integrals
            res = [(i, r) for i, r in enumerate(self.results, 1) if r.channels_a]
            if res:
                p = folder / f"{prefix}_orbital_integrals.csv"
                with open(p, "w", newline="", encoding="utf-8-sig") as fh:
                    w = csv.writer(fh)
                    w.writerow(["No", "bond", "orbital_A", "orbital_B", "quantity", "spin", "integral"])
                    for i, r in res:
                        for q in qs:
                            for sp in spins:
                                M = r.resolved_integrals[q][sp]
                                for ia, (na, _) in enumerate(r.channels_a):
                                    for ib, (nb, _) in enumerate(r.channels_b):
                                        w.writerow([i, r.bond.label, na, nb, QLABEL[q], spin_names[sp], f"{M[ia, ib]:.6f}"])
                written.append(p)
            if lobster_files:
                ok = [r for r in self.results if r.status not in ("aliased", "nonperiodic") and r.bond.mode == "IMAGE"]
                bonds = [r.bond for r in ok]
                for q in qs:
                    kind = QLABEL[q]
                    curves, cumul = [], []
                    for r in ok:
                        _, c = self.curves(r, q, grid=g)
                        _, ci = self.curves(r, q, cumulative=True, grid=g)
                        # restricted: LOBSTER files contain the total incl. both spins already
                        curves.append(c)
                        cumul.append(ci)
                    p = folder / f"{kind}CAR.lobster"
                    lobster.write_coxxcar(p, kind, g, self.eref, bonds, curves, cumul, spins, st)
                    written.append(p)
                    p = folder / f"I{kind}LIST.lobster"
                    integ = [{sp: self._spin_integral(r, q, sp) for sp in spins} for r in ok]
                    lobster.write_icoxxlist(p, kind, bonds, integ, spins, st, sysm.ns)
                    written.append(p)
        if lobster_files:
            p = folder / "POSCAR"
            lobster.write_poscar(p, st)
            written.append(p)
        # populations
        if sysm.has_S:
            if self.pop_m is None:
                self.compute_populations(progress)
            p = folder / f"{prefix}_charges.csv"
            with open(p, "w", newline="", encoding="utf-8-sig") as fh:
                w = csv.writer(fh)
                head = ["atom", "element", "Z_eff", "Mulliken_pop", "Mulliken_charge", "Loewdin_pop", "Loewdin_charge"]
                if sysm.ns == 2:
                    head += ["Mulliken_moment", "Loewdin_moment"]
                w.writerow(head)
                for i, a in enumerate(sysm.d.atoms):
                    row = [a.index, a.element, a.z_eff, f"{self.pop_m.atom_pop[:, i].sum():.5f}",
                           f"{self.pop_m.charges[i]:.5f}", f"{self.pop_l.atom_pop[:, i].sum():.5f}",
                           f"{self.pop_l.charges[i]:.5f}"]
                    if sysm.ns == 2:
                        row += [f"{self.pop_m.moments[i]:.5f}", f"{self.pop_l.moments[i]:.5f}"]
                    w.writerow(row)
            written.append(p)
            p = folder / f"{prefix}_grosspop.csv"
            with open(p, "w", newline="", encoding="utf-8-sig") as fh:
                w = csv.writer(fh)
                w.writerow(["atom", "element", "shell", "l", "orbital", "Mulliken_up_or_total", "Mulliken_down", "Mulliken_total"])
                for row in gross_populations(sysm, self.pop_m):
                    w.writerow([*row[:5], f"{row[5]:.6f}", "" if row[6] is None else f"{row[6]:.6f}", f"{row[7]:.6f}"])
            written.append(p)
            if lobster_files:
                p = folder / "CHARGE.lobster"
                lobster.write_charge(p, st, self.pop_m.charges, self.pop_l.charges)
                written.append(p)
                p = folder / "GROSSPOP.lobster"
                lobster.write_grosspop(p, sysm, self.pop_m, self.pop_l)
                written.append(p)
            if pdos:
                groups = [PDOSGroup(f"{a.element}{a.index}", [a.index], "lm") for a in sysm.d.atoms]
                _, cur = self.pdos(groups, "mulliken", True, progress=progress)
                p = folder / f"{prefix}_atom_pdos.csv"
                names = list(cur[0].keys())
                with open(p, "w", newline="", encoding="utf-8-sig") as fh:
                    w = csv.writer(fh)
                    head = ["E-Eref_eV"] + [f"{n}[{spin_names[s]}]" for s in cur for n in names]
                    w.writerow(head)
                    cols = [g] + [cur[s][n] for s in cur for n in names]
                    for row in np.column_stack(cols):
                        w.writerow([f"{x:.6g}" for x in row])
                written.append(p)
                if lobster_files:
                    pc = {}
                    for a in sysm.d.atoms:
                        cols = []
                        for name, idx in lobster._lobster_orbital_columns(sysm, a.index):
                            cols.append((name, idx))
                        pc[a.index] = cols
                    tdos = {s: cur[s]["TDOS"] for s in cur}
                    pcurves = self._orbital_curves(pc, g)
                    p = folder / "DOSCAR.lobster"
                    lobster.write_doscar(p, sysm, g, tdos, pcurves, list(cur.keys()))
                    written.append(p)
        # report
        p = folder / f"{prefix}_report.json"
        p.write_text(json.dumps(self.report(written), ensure_ascii=False, indent=2), encoding="utf-8")
        written.append(p)
        return written

    def _spin_integral(self, r, q, sp):
        sysm = self.system
        occ = np.array([sysm.occ(k, sp) for k in range(sysm.nk)])
        return float(np.sum(sysm.wts[:, None] * occ * r.values[q][sp]))

    def _orbital_curves(self, pc, grid):
        """Mulliken PDOS for explicit AO subsets: pc[atom] = [(name, ao_idx)] -> [(name, {s: curve})]."""
        from .population import state_weights
        from .spectra import broaden
        sysm = self.system
        s_ = self.settings
        plan = [(atom, name, idx) for atom, cols in pc.items() for name, idx in cols]
        out = {atom: [] for atom in pc}
        res = {}
        for s in range(sysm.ns):
            E, Wt = [], []
            for k in range(sysm.nk):
                W = state_weights(sysm, k, s, "mulliken")
                fac = sysm.g * sysm.wts[k]
                E.append(sysm.eps_ev(k, s) - self.eref)
                Wt.append(np.column_stack([W[idx].sum(axis=0) * fac for _, _, idx in plan]))
            cur = broaden(np.concatenate(E), np.vstack(Wt), grid, s_.fwhm, s_.shape)
            res[s] = cur
        for j, (atom, name, _) in enumerate(plan):
            out[atom].append((name, {s: res[s][:, j] for s in res}))
        return out

    def report(self, written=()):
        sysm = self.system
        info = self.summary()
        inputs = {}
        for key in ("mokp",):
            pth = self.inputs.get(key)
            if pth and Path(pth).exists():
                h = hashlib.sha256()
                with open(pth, "rb") as fh:
                    for chunk in iter(lambda: fh.read(1 << 20), b""):
                        h.update(chunk)
                inputs[key] = {"path": pth, "sha256": h.hexdigest()}
        return {
            "program": f"CP2K-BondDOS {__version__}",
            "time": _dt.datetime.now().isoformat(timespec="seconds"),
            "inputs": {**self.inputs, **inputs},
            "settings": {**self.settings.__dict__, "sigma_eV": self.settings.fwhm / FWHM_PER_SIGMA},
            "system": info,
            "checks": [{"name": c.name, "value": float(c.value), "ok": bool(c.ok), "note": c.note} for c in sysm.checks],
            "gauge_check": {k: (float(v) if isinstance(v, (float, np.floating)) else (bool(v) if isinstance(v, (bool, np.bool_)) else v))
                            for k, v in self.gauge.items()},
            "energy_reference_eV": self.eref, "energy_reference": self.eref_note,
            "bonds": [{"label": r.bond.label, "a": r.bond.a, "b": r.bond.b, "R_raw": list(r.bond.r_raw),
                       "R_cp2k_wrapped": list(r.bond.r_wrap), "length_A": r.bond.length, "mode": r.bond.mode,
                       "status": r.status, "notes": r.notes,
                       "integrals": {QLABEL[q]: v for q, v in r.integrals.items()}} for r in self.results],
            "outputs": [str(p) for p in written],
            "conventions": {
                "COHP": "raw COHP, negative = bonding; -ICOHP positive = bonding; units eV",
                "COOP": "positive = bonding; units electrons",
                "COBI": "positive = bonding; Löwdin basis = Wiberg-type (LOBSTER-like), AO basis = Mayer-type",
                "integrals": "sum_k w_k f_nk v_nk with CP2K occupations (independent of broadening/window)",
                "spectra": "Gaussian/Lorentzian broadened with FWHM; restricted spectra include spin degeneracy 2",
            },
        }


# ------------------------------------------------------------------ cohp.in

def parse_config(path):
    path = Path(path)
    cfg = {"BOND": []}
    for raw in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        k, v = [x.strip() for x in line.split("=", 1)]
        k = k.upper()
        if k == "BOND":
            cfg["BOND"].append(v)
        else:
            cfg[k] = v
    cfg["_dir"] = str(path.parent)
    return cfg


def _truthy(v, default=True):
    if v is None:
        return default
    return str(v).strip().upper() in ("T", "TRUE", "YES", "Y", "1", ".TRUE.")


def session_from_config(cfg, progress=None):
    base = Path(cfg.get("_dir", "."))

    def rel(p):
        if p is None:
            return None
        q = Path(p)
        return q if q.is_absolute() else base / q

    sess = BondSession()
    s = sess.settings
    s.basis = cfg.get("COHP_BASIS", "AO").lower()
    qs = ["cohp"]
    if _truthy(cfg.get("COOP"), True):
        qs.append("coop")
    if _truthy(cfg.get("COBI"), True):
        qs.append("cobi")
    s.quantities = tuple(qs)
    s.resolution = cfg.get("ORBITAL_RESOLUTION", "L").lower()
    s.spins = cfg.get("SPIN", "AUTO").lower()
    if "FWHM_EV" in cfg:
        s.fwhm = float(cfg["FWHM_EV"])
    elif "SIGMA_EV" in cfg:
        s.fwhm = float(cfg["SIGMA_EV"]) * FWHM_PER_SIGMA
    s.emin = float(cfg.get("EMIN_EV", -15))
    s.emax = float(cfg.get("EMAX_EV", 10))
    s.de = float(cfg.get("DE_EV", 0.01))
    ef = cfg.get("EF_EV", "auto").strip().lower()
    if ef in ("auto", "hoco", "midgap", "fermi"):
        s.eref_mode = ef
    else:
        s.eref_mode = "manual"
        s.eref_value = float(ef)
    ks_dir = rel(cfg.get("KS_DIR")) if cfg.get("KS_DIR") else None
    electrons = float(cfg["ELECTRONS"]) if "ELECTRONS" in cfg else None
    sess.load(rel(cfg["MOKP"]), ks_dir, cfg.get("KS_PATTERN"), cfg.get("CSR_TRIANGLE", "FULL"),
              cfg.get("PERIODIC", "XYZ"), electrons, progress)
    mode = cfg.get("MODE", "IMAGE").upper()
    reqs = [BondRequest.parse(b, mode) for b in cfg["BOND"]]
    auto = None
    if "AUTO_BONDS" in cfg:
        t = cfg["AUTO_BONDS"].split()
        auto = {"dmin": float(t[0]), "dmax": float(t[1])}
        els = [x for x in t[2:]]
        if els:
            auto["elements_a"] = [els[0]]
            auto["elements_b"] = els[1:2] or None
    return sess, reqs, auto


def write_config(path, sess: BondSession, bonds_text: list[str], auto=None, prefix="results/CP2K"):
    s = sess.settings
    inp = sess.inputs
    lines = ["# cohp.in written by CP2K-BondDOS", f"MOKP = {inp.get('mokp')}"]
    if inp.get("ks_dir"):
        lines.append(f"KS_DIR = {inp['ks_dir']}")
    if inp.get("ks_set"):
        lines.append(f"KS_PATTERN = {inp['ks_set']}")
    lines += [f"PERIODIC = {inp.get('periodic', 'XYZ')}", f"CSR_TRIANGLE = {inp.get('triangle', 'FULL')}",
              f"COHP_BASIS = {s.basis.upper()}", f"COOP = {'T' if 'coop' in s.quantities else 'F'}",
              f"COBI = {'T' if 'cobi' in s.quantities else 'F'}", f"ORBITAL_RESOLUTION = {s.resolution.upper()}",
              f"SPIN = {s.spins.upper()}", f"FWHM_EV = {s.fwhm:g}", f"EMIN_EV = {s.emin:g}", f"EMAX_EV = {s.emax:g}",
              f"DE_EV = {s.de:g}",
              f"EF_EV = {s.eref_value if s.eref_mode == 'manual' else s.eref_mode}",
              f"OUTPUT_PREFIX = {prefix}"]
    for b in bonds_text:
        lines.append(f"BOND = {b}")
    if auto:
        lines.append(f"AUTO_BONDS = {auto['dmin']} {auto['dmax']}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
