"""Readers for CP2K .dos / .pdos files (old Gamma format and the 2026.2 &DOS interface).

Supported layouts
-----------------
* raw per-MO PDOS   : ``#  MO Eigenvalue [a.u.] Occupation [Total] s py pz px d-2 ...``
* broadened curves  : ``# Energy[eV] total s p d`` (``*_curve*.pdos``, incl. k-point PDOS)
* total DOS         : 3 columns (E, Density, Occupation) or 5/6 columns for two spins
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .constants import HARTREE_EV, CP2K_COMPONENT_ALIASES, L_LETTER

_NUM = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
_ENERGY_LABEL = re.compile(r"(Eigenvalue|Energy)\s*\[\s*([^\]]+?)\s*\]", re.I)


class DOSFileError(ValueError):
    pass


def _unit_factor(unit: str) -> float:
    key = unit.lower().replace(".", "").replace(" ", "")
    if key in ("au", "ha", "hartree", "hartrees"):
        return HARTREE_EV
    if key == "ev":
        return 1.0
    raise DOSFileError(f"不支持的能量单位 {unit!r}")


def _numbers(text):
    return [float(x.replace("D", "E").replace("d", "e")) for x in re.findall(_NUM, text)]


@dataclass(eq=False)
class DosData:
    """One spectrum (one KIND/list, one spin) read from a CP2K file. Energies are stored in eV."""
    path: str
    kind: str                      # 'TDOS', KIND name or 'list<n>'
    spin: str                      # 'none' | 'alpha' | 'beta'
    mode: str                      # 'sticks' (per MO) | 'curve' (broadened) | 'bins' (raw histogram)
    energy: np.ndarray             # eV; absolute for sticks/bins, shifted by `zero` for curves
    channels: dict                 # name -> values (curve: states/eV; sticks/bins: weights)
    zero: str = "absolute"         # 'absolute' | 'fermi' | 'hoco' | 'unknown'
    fermi: float | None = None     # absolute eV
    hoco: float | None = None      # absolute eV
    occupation: np.ndarray | None = None
    fwhm: float | None = None      # broadening of a curve file (eV)
    kpoint: bool = False
    iteration: int | None = None
    components: list = field(default_factory=list)   # m-resolved channel names present
    header: list = field(default_factory=list)

    @property
    def name(self):
        return Path(self.path).name

    @property
    def label(self):
        sp = "" if self.spin == "none" else f"[{self.spin}]"
        return f"{self.kind}{sp}"

    def reference_value(self, ref: str):
        """Absolute eV value of a reference level, or None if unknown."""
        ref = ref.lower()
        if ref == "absolute":
            return 0.0
        if ref == "fermi":
            return self.fermi
        if ref == "hoco":
            return self.hoco
        raise ValueError(ref)

    def energies_relative_to(self, ref: str):
        """Energies on the requested zero (eV). Returns None when the conversion is impossible."""
        ref = ref.lower()
        if ref == "file":
            return self.energy.copy()
        if self.zero == ref:
            return self.energy.copy()
        cur = 0.0 if self.zero == "absolute" else self.reference_value(self.zero) if self.zero in ("fermi", "hoco") else None
        new = self.reference_value(ref)
        if cur is None or new is None:
            return None
        return self.energy + cur - new

    def l_channels(self):
        """Channels collapsed to l (s/p/d/f) plus total; m-components kept as well."""
        return self.channels


def _split_blocks(lines):
    """APPEND files contain several blocks; a new title line after data starts a block."""
    blocks, hdr, rows = [], [], []
    for ln in lines:
        s = ln.strip()
        if not s:
            continue
        if s.startswith("#"):
            starts = re.search(r"iteration step", s, re.I) is not None
            if rows and starts:
                blocks.append((hdr, rows))
                hdr, rows = [], []
            hdr.append(s.lstrip("#").strip())
        else:
            rows.append(s)
    if rows:
        blocks.append((hdr, rows))
    if not blocks:
        raise DOSFileError("文件中没有数值数据")
    return blocks


def _ref_values(headers, name):
    vals = []
    for h in headers:
        m = re.search(rf"E\(\s*{name}\s*\)\s*=(.*)", h, re.I)
        if m:
            tail = m.group(1)
            u = re.search(r"(a\.u\.|Hartree|eV)", tail, re.I)
            if u:
                vals = [v * _unit_factor(u.group(1)) for v in _numbers(tail[:u.start()])]
            else:
                vals = [v * HARTREE_EV for v in _numbers(tail)]
            break
        if name.lower() == "fermi":
            m = re.search(r"E_Fermi\s*\[([^\]]+)\]\s*=\s*(.*)", h, re.I)
            if m:
                vals = [v * _unit_factor(m.group(1)) for v in _numbers(m.group(2))]
                break
    return vals


def _parse_rows(rows):
    data = []
    ncol = None
    for i, r in enumerate(rows):
        toks = r.split()
        try:
            vals = [float(t.replace("D", "E").replace("d", "e")) for t in toks]
        except ValueError:
            raise DOSFileError(f"第 {i + 1} 行数据无法解析（可能是 CP2K 输出 '*****' 溢出）: {r[:60]}")
        if ncol is None:
            ncol = len(vals)
        elif len(vals) != ncol:
            raise DOSFileError("数据列数不一致，文件可能被截断或混入了多次输出")
        data.append(vals)
    return np.asarray(data, float)


def _collapse_components(comp: dict) -> dict:
    """From m-resolved channels build s/p/d/f sums (keeps the components)."""
    out = {}
    for lch in L_LETTER[:4]:
        keys = [k for k in comp if (k == "s" and lch == "s") or (k != "s" and k[0] == lch)]
        if keys:
            out[lch] = np.sum([comp[k] for k in keys], axis=0)
    return out


def _channels_from_labels(labels, cols):
    """labels: raw CP2K column labels (after energy/occupation)."""
    low = [x.strip().lower() for x in labels]
    chans = {}
    comps = []
    total = None
    for lab, col in zip(low, cols):
        if lab == "total":
            total = col
            continue
        if lab in ("s", "p", "d", "f", "g"):
            chans[lab] = chans.get(lab, 0) + col
            continue
        friendly = CP2K_COMPONENT_ALIASES.get(lab, lab)
        chans[friendly] = col
        comps.append(friendly)
    if comps:
        chans.update({k: v for k, v in _collapse_components({c: chans[c] for c in comps}).items() if k not in chans})
    if total is None:
        ls = [chans[k] for k in ("s", "p", "d", "f", "g") if k in chans]
        if not ls:
            ls = [chans[c] for c in comps]
        total = np.sum(ls, axis=0) if ls else np.zeros_like(cols[0])
    chans = {"total": total, **chans}
    return chans, comps


def read_dos_file(path, block: int | str = "last") -> list[DosData]:
    """Read a CP2K .dos/.pdos file and return one DosData per spin channel."""
    path = Path(path)
    try:
        text = path.read_text(encoding="utf-8-sig", errors="replace")
    except OSError as e:
        raise DOSFileError(f"无法读取 {path}: {e}")
    blocks = _split_blocks(text.splitlines())
    if block == "last":
        hdr, rows = blocks[-1]
    else:
        hdr, rows = blocks[int(block) - 1]
    data = _parse_rows(rows)
    if not np.all(np.isfinite(data)):
        raise DOSFileError("数据中含 NaN/Inf")
    col_line = None
    for h in reversed(hdr):
        if _ENERGY_LABEL.search(h):
            col_line = h
            break
    if col_line is None:
        raise DOSFileError(f"{path.name}: 找不到能量列标题（Energy[...] / Eigenvalue [...]）")
    m = _ENERGY_LABEL.search(col_line)
    factor = _unit_factor(m.group(2))
    after = col_line[m.end():].split()
    before = col_line[:m.start()].split()
    joined = "\n".join(hdr)
    fermi_v = _ref_values(hdr, "Fermi")
    hoco_v = _ref_values(hdr, "HOCO")
    spin_m = re.search(r"(?:^|[-_.])(ALPHA|BETA)(?:[-_.]|$)", path.name, re.I)
    spin = spin_m.group(1).lower() if spin_m else "none"
    zero_line = [h for h in hdr if re.match(r"Energy zero\s*:", h, re.I)]
    zero = zero_line[-1].split(":", 1)[1].split("->")[-1].strip().lower() if zero_line else None
    fw = None
    for h in hdr:
        mm = re.search(r"broadening,\s*FWHM\s*=\s*(" + _NUM + r")\s*a\.u\.\s*=\s*(" + _NUM + r")\s*eV", h, re.I)
        if mm:
            fw = float(mm.group(2))
    it = None
    mm = re.search(r"iteration step\s*i\s*=\s*(\d+)", joined, re.I)
    if mm:
        it = int(mm.group(1))
    kpoint = re.search(r"K-point", joined, re.I) is not None
    km = re.search(r"atomic kind\s+(\S+)", joined, re.I)
    lm = re.search(r"(?:atom list|for list)\s+(\d+)", joined, re.I)
    is_tdos = re.search(r"DOS at iteration", joined, re.I) is not None and km is None and lm is None
    kind = km.group(1) if km else (f"list{lm.group(1)}" if lm else ("TDOS" if is_tdos else path.stem))
    is_raw_mo = any(t.upper() == "MO" for t in before)

    def pick(vals, sp):
        if not vals:
            return None
        if len(vals) >= 2 and sp == "beta":
            return vals[1]
        return vals[0]

    common = dict(path=str(path), kind=kind, header=hdr, iteration=it, kpoint=kpoint)
    if is_raw_mo:
        # columns: MO, E, Occupation, [Total], components...
        labels = after
        if not labels or labels[0].lower() != "occupation":
            raise DOSFileError(f"{path.name}: 原始 PDOS 表头应在能量后为 Occupation")
        labels = labels[1:]
        if len(labels) != data.shape[1] - 3:
            raise DOSFileError(f"{path.name}: 表头列数 ({len(labels)}) 与数据列数 ({data.shape[1] - 3}) 不符")
        chans, comps = _channels_from_labels(labels, [data[:, 3 + i] for i in range(len(labels))])
        e_abs = data[:, 1] * factor
        occ = data[:, 2]
        hoco = pick(hoco_v, spin)
        if hoco is None and np.any(occ > 1e-8):
            hoco = float(e_abs[occ > 1e-8].max())
        return [DosData(spin=spin, mode="sticks", energy=e_abs, channels=chans, zero="absolute",
                        fermi=pick(fermi_v, spin), hoco=hoco, occupation=occ, components=comps, **common)]
    if is_tdos:
        curve = zero is not None or fw is not None
        if not curve:
            zero = "absolute"
        layouts = {3: [(0, 1, 2, spin)], 5: [(0, 1, 2, "alpha"), (0, 3, 4, "beta")],
                   6: [(0, 1, 2, "alpha"), (3, 4, 5, "beta")]}
        if data.shape[1] not in layouts:
            raise DOSFileError(f"{path.name}: total DOS 应为 3/5/6 列，实际 {data.shape[1]} 列")
        out = []
        for ie, idn, io, sp in layouts[data.shape[1]]:
            mask = np.ones(len(data), bool) if curve else (data[:, idn] != 0) | (data[:, io] != 0)
            dens = data[mask, idn] / (factor if curve else 1.0)
            occ = data[mask, io] / (factor if curve else 1.0)
            out.append(DosData(spin=sp, mode="curve" if curve else "bins", energy=data[mask, ie] * factor,
                               channels={"total": dens, "occupied": occ}, zero=zero or "unknown",
                               fermi=pick(fermi_v, sp), hoco=pick(hoco_v, sp), fwhm=fw, **common))
        return out
    # broadened PDOS curve
    labels = after
    if len(labels) != data.shape[1] - 1:
        raise DOSFileError(f"{path.name}: 曲线表头列数 ({len(labels)}) 与数据列数 ({data.shape[1] - 1}) 不符")
    cols = [data[:, 1 + i] / factor for i in range(len(labels))]
    chans, comps = _channels_from_labels(labels, cols)
    return [DosData(spin=spin, mode="curve", energy=data[:, 0] * factor, channels=chans,
                    zero=zero or "unknown", fermi=pick(fermi_v, spin), hoco=pick(hoco_v, spin), fwhm=fw,
                    components=comps, **common)]


def scan_folder(folder):
    """Classify CP2K output files found in a folder (non recursive)."""
    folder = Path(folder)
    found = {"dos": [], "pdos": [], "mokp": [], "csr": [], "out": []}
    if not folder.is_dir():
        return found
    for p in sorted(folder.iterdir()):
        if not p.is_file():
            continue
        n = p.name.lower()
        if n.endswith(".pdos"):
            found["pdos"].append(p)
        elif n.endswith(".dos"):
            found["dos"].append(p)
        elif n.endswith(".mokp"):
            found["mokp"].append(p)
        elif n.endswith(".csr"):
            found["csr"].append(p)
        elif n.endswith((".out", ".log")):
            found["out"].append(p)
    return found
