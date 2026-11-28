"""COHP / COOP / COBI engine for CP2K k-point data.

Conventions
-----------
H(k) = sum_R H(R) exp(+2πi k·R)  (CP2K rskp_transform, R = integer cell of the neighbour atom in the
wrapped frame used by CP2K's neighbour lists).  For a bond A(0)–B(R):

    M̃_AB(R)   = Σ_k w_k exp(-2πi k·R) M_AB(k)
    v_nkσ     = 2 Re[ x_A,nk† M̃_AB(R) exp(+2πi k·R) y_B,nk ]          (IMAGE mode)
    v_nkσ     = 2 Re[ x_A,nk† M_AB(k) y_B,nk ]                          (PAIR mode, all images summed)

    quantity  x      y     M
    COHP-AO   C      C     H                (direct Gaussian AO)
    COHP-L    Q      Q     S^-1/2 H S^-1/2  (Löwdin, Q = S^1/2 C)
    COOP      C      C     S
    COBI-L    Q      Q     P^σ_L = Q diag(f/g) Q†            (Wiberg type, LOBSTER-like)
    COBI-AO   S C    C     (P^σ S) = C diag(f/g) (SC)†      (Mayer type)

Gauge of the AO-basis COHP (``cohp_gauge``): in a non-orthogonal basis C†HC changes by Δ·COOP under a
uniform shift H → H + ΔS of the (arbitrary) potential zero of a periodic calculation.
    ecov   H_μν − ½(ε_μ + ε_ν) S_μν,  ε_μ = H_μμ(R=0)   (DEFAULT; covalent bond energy, shift invariant;
           Börnsen, Meyer, Grotheer, Fähnle, J. Phys.: Condens. Matter 11, L287 (1999))
    eref   H_μν − E_ref S_μν                      (energies measured from the chosen energy zero)
    raw    H_μν                                   (v1.0.0 behaviour; depends on the potential zero)
The Löwdin COHP is shift invariant by construction and ignores this setting.

Spectrum:  X(E) = Σ g w_k v_nk δ(E − ε_nk);   integral  IX = Σ w_k f_nk v_nk  (f = CP2K occupation).
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np

from .constants import l_name, m_name
from .kspace import KSystem, Bond
from .spectra import broaden, make_grid

QUANTITIES = ("cohp", "coop", "cobi")
QLABEL = {"cohp": "COHP", "coop": "COOP", "cobi": "COBI"}
QUNIT = {"cohp": "eV", "coop": "e", "cobi": ""}


# ------------------------------------------------------------------ orbital channels

def atom_channels(system: KSystem, atom: int, resolution: str):
    """List of (name, local_indices) for an atom. resolution: none|l|lm|shell|ao."""
    d = system.d
    idx = system.atom_aos(atom)
    ls = d.ao_l[idx]
    ms = d.ao_m[idx]
    sh = d.ao_shell[idx]
    res = (resolution or "none").lower()
    if res == "none":
        return [("all", np.arange(idx.size))]
    groups = {}
    order = []
    for j in range(idx.size):
        if res == "l":
            key = l_name(ls[j])
        elif res == "lm":
            key = m_name(ls[j], ms[j])
        elif res == "shell":
            key = f"{l_name(ls[j])}({sh[j]})"
        else:  # ao
            key = f"{m_name(ls[j], ms[j])}({sh[j]})"
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append(j)
    return [(k, np.array(groups[k])) for k in order]


def indicator(channels, n):
    g = np.zeros((len(channels), n))
    for i, (_, loc) in enumerate(channels):
        g[i, loc] = 1.0
    return g


# ------------------------------------------------------------------ results

@dataclass(eq=False)
class BondResult:
    bond: Bond
    status: str = "ok"
    notes: list = field(default_factory=list)
    # values[q][s] -> (nk, nmo) per-state values (total)
    values: dict = field(default_factory=dict)
    # resolved[q][s] -> (nk, nchA, nchB, nmo)
    resolved: dict = field(default_factory=dict)
    channels_a: list = field(default_factory=list)
    channels_b: list = field(default_factory=list)
    integrals: dict = field(default_factory=dict)        # q -> {'up','down','total'} -> float
    resolved_integrals: dict = field(default_factory=dict)  # q -> s -> (nchA, nchB)

    @property
    def label(self):
        return self.bond.label


@dataclass
class BondingOptions:
    basis: str = "ao"                      # ao | lowdin
    quantities: tuple = ("cohp", "coop", "cobi")
    resolution: str = "none"               # none|l|lm|shell|ao  (applied to bonds marked resolve=True)
    spins: str = "auto"                    # auto|up|down|both
    cohp_gauge: str = "ecov"               # ecov|eref|raw  (AO basis only)
    cohp_shift: float = 0.0                # eV, reference energy used by cohp_gauge = "eref"


class BondingEngine:
    def __init__(self, system: KSystem, options: BondingOptions):
        self.sys = system
        self.opt = options
        qs = [q for q in options.quantities if q in QUANTITIES]
        if "cohp" in qs and not system.has_H:
            qs.remove("cohp")
        if not system.has_S:
            qs = [q for q in qs if q == "cohp"]
            if options.basis == "lowdin":
                raise ValueError("Löwdin 表象需要重叠矩阵 S(k)（MO_KP 需 AO_EXPORT_TYPE OVERLAP_MATRIX）")
        if not qs:
            raise ValueError("没有可计算的量：COHP 需要 KS 矩阵，COOP/COBI 需要 S(k)")
        self.quantities = qs
        ns = system.ns
        sp = (options.spins or "auto").lower()
        if ns == 1:
            self.spins = [0]
        elif sp == "up":
            self.spins = [0]
        elif sp == "down":
            self.spins = [1]
        else:
            self.spins = [0, 1]
        self.log = []
        self.gauge = (options.cohp_gauge or "ecov").lower()
        self.gauge_fallback = False
        if self.gauge not in ("raw", "ecov", "eref"):
            raise ValueError(f"未知的 COHP_GAUGE: {options.cohp_gauge!r}（可选 RAW | ECOV | EREF）")
        if options.basis == "lowdin" or "cohp" not in qs:
            self.gauge = "raw"
        elif self.gauge != "raw" and not system.has_S:
            self.gauge, self.gauge_fallback = "raw", True      # no S(k): only the raw gauge is possible
        self._onsite_cache = {}

    def onsite_energies(self, s):
        """ε_μ = H_μμ(R=0) = Σ_k w_k H_μμ(k) in eV for spin s (used by the E_cov gauge)."""
        if s not in self._onsite_cache:
            sysm = self.sys
            self._onsite_cache[s] = sum(sysm.wts[k] * np.real(np.diag(sysm.H_ev(k, s))) for k in range(sysm.nk))
        return self._onsite_cache[s]

    # ---- bond construction
    def make_bond(self, a, b, r_raw=None, mode="IMAGE", all_images=False):
        st = self.sys.structure
        n = len(st.atoms)
        if not (1 <= a <= n and 1 <= b <= n):
            raise ValueError(f"原子编号超出 1..{n}: {a} {b}")
        mode = mode.upper()
        if mode == "PAIR":
            if a == b:
                raise ValueError("PAIR 模式需要两个不同原子")
            R, d = st.nearest(a, b)[0]
            return [Bond(a, b, (0, 0, 0), (0, 0, 0), d, "PAIR", f"{st.label(a)}–{st.label(b)} (all images)")]
        if r_raw is None:
            sel = st.nearest(a, b, all_equivalent=all_images)
        else:
            r_raw = tuple(int(x) for x in r_raw)
            if a == b and r_raw == (0, 0, 0):
                raise ValueError("同一原子 R=(0,0,0) 不是化学键")
            sel = [(r_raw, st.distance(a, b, r_raw))]
        out = []
        for R, d in sel:
            rw = st.wrap_translation(a, b, R)
            out.append(Bond(a, b, tuple(R), rw, d, "IMAGE", f"{st.label(a)}–{st.label(b)} [{R[0]} {R[1]} {R[2]}]"))
        return out

    # ---- per-k data
    def _k_payload(self, k, s, atoms, qs):
        """Return {q: (XL, XR, Mfull_or_None, Mfn)} restricted later by atom blocks."""
        sysm = self.sys
        C = sysm.C(k, s)
        f = sysm.occ(k, s) / sysm.g
        out = {}
        S = sysm.S(k) if sysm.has_S else None
        need_lowdin = self.opt.basis == "lowdin" and any(q in qs for q in ("cohp", "cobi"))
        if need_lowdin:
            sh, sih, _, _ = sysm.lowdin(k)
            Q = sh @ C
        if "cohp" in qs:
            H = sysm.H_ev(k, s)
            if self.opt.basis == "lowdin":
                out["cohp"] = (Q, Q, sih @ H @ sih)
            else:
                if self.gauge == "ecov":
                    e = self.onsite_energies(s)
                    H = H - 0.5 * (e[:, None] * S + S * e[None, :])
                elif self.gauge == "eref":
                    H = H - self.opt.cohp_shift * S
                out["cohp"] = (C, C, H)
        if "coop" in qs:
            out["coop"] = (C, C, S)
        if "cobi" in qs:
            if self.opt.basis == "lowdin":
                out["cobi"] = (Q, Q, ("low", Q, f))
            else:
                D = S @ C
                out["cobi"] = (D, C, ("mayer", C, D, f))
        return out

    @staticmethod
    def _block(M, ia, ib):
        if isinstance(M, tuple):
            if M[0] == "low":
                _, Q, f = M
                return (Q[ia] * f[None, :]) @ Q[ib].conj().T
            _, C, D, f = M
            return (C[ia] * f[None, :]) @ D[ib].conj().T
        return M[np.ix_(ia, ib)]

    def run(self, bonds: list[Bond], resolve: set | None = None, progress=None, cancel=None):
        """Compute per-state values for all bonds. resolve: set of bond indices to orbital-resolve."""
        t0 = time.time()
        sysm = self.sys
        qs = self.quantities
        resolve = resolve or set()
        results = [BondResult(b) for b in bonds]
        aos = {}
        for b in bonds:
            for at in (b.a, b.b):
                if at not in aos:
                    aos[at] = sysm.atom_aos(at)
        for i, r in enumerate(results):
            b = r.bond
            if b.mode == "IMAGE":
                st, notes = sysm.image_status(b.r_wrap)
                if not sysm.mesh.full and sysm.nk > 1:
                    st, notes = "aliased", ["k 点网格不是完整等权网格，无法做周期镜像分解"]
                r.status, r.notes = st, notes
            elif not sysm.mesh.full and sysm.nk > 1:
                r.notes = ["k 网格经过点群对称性约化：PAIR 结果只有对对称等价的键取平均后才有意义（请用 SYMMETRY F）"]
            if i in resolve:
                r.channels_a = atom_channels(sysm, b.a, self.opt.resolution)
                r.channels_b = atom_channels(sysm, b.b, self.opt.resolution)
        kp = sysm.kpts
        w = sysm.wts
        image_idx = [i for i, r in enumerate(results) if r.bond.mode == "IMAGE" and r.status not in ("aliased", "nonperiodic")]
        # pass 1: Fourier transform M blocks for IMAGE bonds
        mt = {}
        total_steps = (2 if image_idx else 1) * sysm.nk * len(self.spins)
        step = 0
        if image_idx:
            for k in range(sysm.nk):
                for s in self.spins:
                    if cancel and cancel():
                        raise RuntimeError("已取消")
                    step += 1
                    if progress:
                        progress(step / total_steps, f"Fourier 变换 H(k)→H(R): k {k + 1}/{sysm.nk}")
                    pay = self._k_payload(k, s, aos, qs)
                    for i in image_idx:
                        b = results[i].bond
                        ph = np.exp(-2j * np.pi * np.dot(kp[k], b.r_wrap)) * w[k]
                        for q in qs:
                            blk = self._block(pay[q][2], aos[b.a], aos[b.b])
                            key = (i, q, s)
                            mt[key] = mt.get(key, 0) + ph * blk
        # pass 2: state contributions
        for r in results:
            for q in qs:
                r.values[q] = {s: np.zeros((sysm.nk, sysm.nmo)) for s in self.spins}
                if r.channels_a:
                    r.resolved[q] = {s: np.zeros((sysm.nk, len(r.channels_a), len(r.channels_b), sysm.nmo))
                                     for s in self.spins}
        GA = {}
        for i, r in enumerate(results):
            if r.channels_a:
                GA[i] = (indicator(r.channels_a, aos[r.bond.a].size), indicator(r.channels_b, aos[r.bond.b].size))
        for k in range(sysm.nk):
            for s in self.spins:
                if cancel and cancel():
                    raise RuntimeError("已取消")
                step += 1
                if progress:
                    progress(step / total_steps, f"计算成键贡献: k {k + 1}/{sysm.nk}")
                pay = self._k_payload(k, s, aos, qs)
                for i, r in enumerate(results):
                    b = r.bond
                    ia, ib = aos[b.a], aos[b.b]
                    if r.status in ("aliased", "nonperiodic"):
                        continue
                    for q in qs:
                        XL, XR, M = pay[q]
                        if b.mode == "IMAGE":
                            Mt = mt[(i, q, s)] * np.exp(2j * np.pi * np.dot(kp[k], b.r_wrap))
                        else:
                            Mt = self._block(M, ia, ib)
                        xa = XL[ia]
                        yb = XR[ib]
                        if i in GA:
                            # (nA, nB, nmo) orbital pair contributions
                            W = 2.0 * np.real(xa.conj()[:, None, :] * Mt[:, :, None] * yb[None, :, :])
                            ga, gb = GA[i]
                            Rk = np.einsum("am,mnk,bn->abk", ga, W, gb, optimize=True)
                            r.resolved[q][s][k] = Rk
                            r.values[q][s][k] = Rk.sum(axis=(0, 1))
                        else:
                            r.values[q][s][k] = 2.0 * np.real(np.sum(xa.conj() * (Mt @ yb), axis=0))
        # integrals
        for r in results:
            for q in qs:
                tot = {}
                for s in self.spins:
                    occ = np.array([sysm.occ(k, s) for k in range(sysm.nk)])
                    tot[s] = float(np.sum(w[:, None] * occ * r.values[q][s]))
                    if r.channels_a:
                        r.resolved_integrals.setdefault(q, {})[s] = np.einsum(
                            "k,kn,kabn->ab", w, occ, r.resolved[q][s])
                r.integrals[q] = _spin_dict(tot, sysm.ns)
        self.log.append(f"成键分析完成：{len(bonds)} 根键，{len(qs)} 种量，用时 {time.time() - t0:.1f} s")
        return results

    # ---- spectra
    def spectra(self, result: BondResult, q, eref, grid, fwhm, shape="gaussian", channel=None, cumulative=False):
        """Return {spin: curve} on grid; channel=(ia, ib) selects an orbital pair (None = total)."""
        sysm = self.sys
        out = {}
        for s in self.spins:
            e = np.concatenate([sysm.eps_ev(k, s) - eref for k in range(sysm.nk)])
            if channel is None:
                v = result.values[q][s]
            else:
                v = result.resolved[q][s][:, channel[0], channel[1], :]
            wt = (sysm.g * sysm.wts[:, None] * v).ravel()
            out[s] = broaden(e, wt, grid, fwhm, shape, cumulative=cumulative)
        return out


def _spin_dict(tot, ns):
    if ns == 1:
        return {"total": tot[0]}
    d = {}
    if 0 in tot:
        d["up"] = tot[0]
    if 1 in tot:
        d["down"] = tot[1]
    d["total"] = sum(tot.values())
    return d


# ------------------------------------------------------------------ gauge diagnostics

def overlap_images(system: KSystem, a, b, rmax=None):
    """S̃_AB(R_wrap) Frobenius norms for all R inside the alias-free box."""
    mesh = system.mesh.mesh
    per = system.structure.periodic
    rng = []
    for ax in range(3):
        n = mesh[ax]
        if not per[ax] or n == 1:
            rng.append([0])
        else:
            h = (n - 1) // 2
            if rmax is not None:
                h = min(h, rmax)
            rng.append(list(range(-h, h + 1)))
    ia, ib = system.atom_aos(a), system.atom_aos(b)
    blocks = [system.S(k)[np.ix_(ia, ib)] for k in range(system.nk)]
    out = []
    for R in np.array(np.meshgrid(*rng, indexing="ij")).reshape(3, -1).T:
        acc = 0
        for k in range(system.nk):
            acc = acc + system.wts[k] * np.exp(-2j * np.pi * np.dot(system.kpts[k], R)) * blocks[k]
        out.append((tuple(int(x) for x in R), float(np.linalg.norm(acc))))
    return out


def gauge_check(system: KSystem, max_atoms=12):
    """Check the Fourier sign convention and wrapped-frame translations with overlap decay.

    For the correct convention the large |S̃_AB(R)| values belong to geometrically short A–B(R) bonds.
    Returns dict with scores; score(+) < score(−) means the implemented convention is consistent.
    """
    st = system.structure
    if not system.has_S or system.nk < 2 or not system.mesh.full:
        return {"available": False, "reason": "需要多 k 点完整网格和 S(k)"}
    n = len(st.atoms)
    atoms = list(range(1, min(n, max_atoms) + 1))
    plus = minus = 0.0
    used = 0
    for a in atoms:
        for b in range(1, n + 1):
            if a == b:
                continue
            for Rw, val in overlap_images(system, a, b, rmax=2):
                if val < 1e-6:
                    continue
                rr = st.raw_translation(a, b, Rw)
                rm = st.raw_translation(a, b, tuple(-x for x in Rw))
                plus += val ** 2 * st.distance(a, b, rr) ** 2
                minus += val ** 2 * st.distance(a, b, rm) ** 2
                used += 1
    if used == 0:
        return {"available": False, "reason": "k 网格太小（每个方向需 ≥3 个点才能区分 ±R）"}
    ratio = minus / plus if plus > 0 else float("inf")
    return {"available": True, "score_plus": plus, "score_minus": minus, "ratio": ratio,
            "consistent": ratio >= 1.0 - 1e-9, "decisive": abs(np.log(max(ratio, 1e-300))) > np.log(1.5)}


def resolve_wrap_ambiguity(system: KSystem):
    """For atoms sitting on a ±0.5 boundary choose the wrap that makes large overlaps short."""
    st = system.structure
    notes = []
    if not st.ambiguous or system.nk < 2 or not system.mesh.full or not system.has_S:
        return notes
    n = len(st.atoms)
    for (atom, ax) in st.ambiguous:
        if system.mesh.mesh[ax] < 3:
            continue
        cands = []
        base = st.g[atom - 1, ax]
        for delta in (0, -1, 1):
            st.g[atom - 1, ax] = base + delta
            score = 0.0
            for b in range(1, n + 1):
                if b == atom:
                    continue
                for Rw, val in overlap_images(system, atom, b, rmax=2):
                    if val < 1e-6:
                        continue
                    score += val ** 2 * st.distance(atom, b, st.raw_translation(atom, b, Rw)) ** 2
            cands.append((score, delta))
        cands.sort()
        st.g[atom - 1, ax] = base + cands[0][1]
        st.frac_wrap[atom - 1] = st.frac_raw[atom - 1] + st.g[atom - 1]
        notes.append(f"{st.label(atom)} 的 {'abc'[ax]} 方向回卷由重叠衰减判据确定 (Δ={cands[0][1]})")
    return notes


def default_grid(system: KSystem, eref, emin=None, emax=None, de=0.01):
    allE = np.concatenate([system.eps_ev(k, s) for (k, s) in system.d.eig]) - eref
    lo = float(allE.min()) - 1.0 if emin is None else emin
    hi = float(allE.max()) + 1.0 if emax is None else emax
    return make_grid(lo, hi, de)
