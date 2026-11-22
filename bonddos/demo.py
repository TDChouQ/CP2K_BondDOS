"""Synthetic non-orthogonal tight-binding model written in the CP2K 2026.2 MO_KP / KS_CSR text formats.

The model is NOT a CP2K calculation; it is used for self tests and for trying the GUI. Its real-space
H(R), S(R) are known exactly so k-space results can be checked against a Γ-point supercell.
"""
from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np

from .constants import HARTREE_EV

# basis per element: list of (shell_index, l)
BASIS = {"C": [(1, 0), (2, 1)], "O": [(1, 0), (2, 0), (3, 1)]}
ONSITE = {"C": {0: -0.55, 1: -0.20}, "O": {0: -0.95, 1: -0.45}}  # Hartree per l (shell 2 of O offset)
ZEFF = {"C": 4, "O": 6}
ZNUC = {"C": 6, "O": 8}


class Model:
    def __init__(self, uks=False, metal=False, seed=20261122, raw_shift=True, cutoff=3.3, gamma_real=False):
        self.cell = np.array([[3.10, 0.00, 0.00], [0.25, 3.00, 0.00], [0.10, -0.15, 3.20]])
        self.elements = ["C", "O", "C"]
        frac_wrapped = np.array([[0.02, 0.03, -0.04], [0.31, 0.27, 0.24], [-0.33, 0.45, 0.12]])
        self.frac_wrapped = frac_wrapped
        # raw (printed) coordinates differ from wrapped ones by lattice vectors -> tests the gauge handling
        shift = np.array([[0, 0, 0], [1, 0, -1], [0, -1, 2]]) if raw_shift else np.zeros((3, 3), int)
        self.frac_raw = frac_wrapped + shift
        self.uks = uks
        self.metal = metal
        self.rng = np.random.default_rng(seed)
        self.cutoff = cutoff
        self.aos = []   # (atom(1-based), shell, l, m)
        for ia, el in enumerate(self.elements, 1):
            for sh, l in BASIS[el]:
                for m in range(-l, l + 1):
                    self.aos.append((ia, sh, l, m))
        self.nao = len(self.aos)
        self._build_real_space()

    def atom_slices(self):
        out = []
        for ia in range(1, len(self.elements) + 1):
            idx = [i for i, a in enumerate(self.aos) if a[0] == ia]
            out.append((idx[0], idx[-1] + 1))
        return out

    def _build_real_space(self):
        """Dict R(tuple) -> (S_R, H_R_list per spin) in the wrapped frame, Hermitian partner included."""
        nao = self.nao
        self.S_R = {}
        self.H_R = [{}, {}]
        sl = self.atom_slices()
        rng = self.rng
        pairs_done = set()
        Rs = list(itertools.product(range(-2, 3), repeat=3))
        for R in Rs:
            self.S_R.setdefault(R, np.zeros((nao, nao)))
            for sp in range(2):
                self.H_R[sp].setdefault(R, np.zeros((nao, nao)))
        for ia in range(3):
            for ib in range(3):
                for R in Rs:
                    key = (ia, ib, R)
                    mkey = (ib, ia, tuple(-x for x in R))
                    if key in pairs_done or mkey in pairs_done:
                        continue
                    pairs_done.add(key)
                    v = (self.frac_wrapped[ib] + np.array(R) - self.frac_wrapped[ia]) @ self.cell
                    d = np.linalg.norm(v)
                    a0, a1 = sl[ia]
                    b0, b1 = sl[ib]
                    if ia == ib and R == (0, 0, 0):
                        blkS = np.eye(a1 - a0)
                        el = self.elements[ia]
                        e = np.array([ONSITE[el][self.aos[j][2]] + (0.12 if self.aos[j][1] == 2 and el == "O" else 0.0)
                                      for j in range(a0, a1)])
                        for sp in range(2):
                            split = (0.03 if sp == 0 else -0.03) if self.uks else 0.0
                            self.H_R[sp][R][a0:a1, a0:a1] += np.diag(e - split * (el == "C"))
                        self.S_R[R][a0:a1, a0:a1] += blkS
                        continue
                    if d > self.cutoff:
                        continue
                    shape = (a1 - a0, b1 - b0)
                    rs = rng.uniform(0.4, 1.0, shape) * rng.choice([-1, 1], shape)
                    rh = rng.uniform(0.5, 1.0, shape)
                    s_blk = 0.10 * np.exp(-1.1 * (d - 1.5)) * rs
                    h_blk = -0.08 * np.exp(-0.9 * (d - 1.5)) * rh * np.sign(rs) + 0.3 * s_blk * (-0.5)
                    self.S_R[R][a0:a1, b0:b1] += s_blk
                    mR = tuple(-x for x in R)
                    self.S_R[mR][b0:b1, a0:a1] += s_blk.T
                    for sp in range(2):
                        hb = h_blk * (1.0 + (0.05 if (self.uks and sp == 1) else 0.0))
                        self.H_R[sp][R][a0:a1, b0:b1] += hb
                        self.H_R[sp][mR][b0:b1, a0:a1] += hb.T

    def hk(self, k, sp):
        nao = self.nao
        H = np.zeros((nao, nao), complex)
        S = np.zeros((nao, nao), complex)
        for R, sR in self.S_R.items():
            ph = np.exp(2j * np.pi * np.dot(k, R))
            S += sR * ph
            H += self.H_R[sp][R] * ph
        return H, S

    def nelectrons(self):
        return sum(ZEFF[e] for e in self.elements)


def mp_mesh(n, gamma_centered=True):
    pts = []
    for i in itertools.product(*[range(x) for x in n]):
        if gamma_centered:
            k = [(i[a] / n[a]) for a in range(3)]
            k = [x - 1.0 if x > 0.5 + 1e-12 else x for x in k]
        else:
            k = [(2 * i[a] - n[a] + 1) / (2 * n[a]) for a in range(3)]
        pts.append(k)
    return np.array(pts)


def solve(model: Model, kpts, kT=0.0):
    """Generalized eigenproblem on all k; returns eig, C, occupations (per spin channel)."""
    def eigh(H, S):
        L = np.linalg.cholesky(S)
        Li = np.linalg.inv(L)
        e, U = np.linalg.eigh(Li @ H @ Li.conj().T)
        return e, Li.conj().T @ U
    ns = 2 if model.uks else 1
    g = 1.0 if ns == 2 else 2.0
    eig, coef = {}, {}
    for ik, k in enumerate(kpts):
        for sp in range(ns):
            H, S = model.hk(k, sp)
            e, C = eigh(H, S)
            eig[(ik, sp)] = e
            coef[(ik, sp)] = C
    nk = len(kpts)
    nel = model.nelectrons()
    w = 1.0 / nk
    allE = np.sort(np.concatenate([eig[key] for key in eig]))
    if kT <= 0:
        # insulating filling: fill the lowest nel/g * nk states (per spin channels pooled)
        nfill = int(round(nel / g * nk)) if ns == 1 else int(round(nel * nk))
        ef = 0.5 * (allE[nfill - 1] + allE[nfill])
        occ = {key: g * (e < ef).astype(float) for key, e in eig.items()}
    else:
        lo, hi = allE[0] - 1, allE[-1] + 1
        for _ in range(200):
            mu = 0.5 * (lo + hi)
            n = sum((w * g / (1 + np.exp(np.clip((e - mu) / kT, -200, 200)))).sum() for e in eig.values())
            lo, hi = (mu, hi) if n < nel else (lo, mu)
        ef = mu
        occ = {key: g / (1 + np.exp(np.clip((e - ef) / kT, -200, 200))) for key, e in eig.items()}
    return eig, coef, occ, ef


def write_mokp(path, model: Model, kpts, eig, coef, occ, ndigits=14, real=False):
    nk = len(kpts)
    ns = 2 if model.uks else 1
    nao = model.nao
    nmo = nao
    w = 1.0 / nk
    sl = model.atom_slices()
    lines = ["# CP2K_KPOINT_MO_DUMP, Version 2.0",
             "# All energies in Hartree, lengths in Angstrom, k-points in fractional reciprocal coords",
             f"# DIMENSIONS: natom nspins nao nkp ={len(model.elements):8d}{ns:8d}{nao:8d}{nk:8d}",
             f"# NMO ={nmo:8d}",
             f"# USE_REAL_WFN = {'T' if real else 'F'}",
             "# AO_EXPORT_TYPE = OVERLAP_MATRIX",
             f"# SPARSE_THRESHOLD ={10.0 ** -ndigits:12.5E}",
             "# CELL_VECTORS [Angstrom]"]
    for row in model.cell:
        lines.append("".join(f"{x:12.6f}   " for x in row))
    lines.append("# ATOM_LIST: Atom_ID  Element  Z_nuc  Z_eff  X [Ang]  Y [Ang]  Z [Ang]  First_AO  Last_AO")
    for ia, el in enumerate(model.elements):
        r = model.frac_raw[ia] @ model.cell
        lines.append(f"{ia + 1:6d}  {el:2s}{ZNUC[el]:6d}{ZEFF[el]:6d}{r[0]:15.6f}{r[1]:15.6f}{r[2]:15.6f}"
                     f"{sl[ia][0] + 1:10d}{sl[ia][1]:10d}")
    lines.append("# KPOINT_LIST: ikp  kx  ky  kz  weight")
    for ik, k in enumerate(kpts):
        lines.append(f"{ik + 1:6d}" + "".join(f"{x:18.8E}" for x in (*k, w)))
    lines.append("# AO_LIST: ao_index  atom_id  atom_shell_index  l  m")
    for i, (ia, sh, l, m) in enumerate(model.aos):
        lines.append(f"{i + 1:8d}{ia:8d}{sh:8d}{l:8d}{m:8d}")
    fw = ndigits + 10
    thresh = 10.0 ** -ndigits

    def vec(v):
        out = []
        for i in range(0, len(v), 5):
            out.append("".join(f"{x:{fw}.{ndigits}E}" for x in v[i:i + 5]))
        return out

    def sparse(M, rowfirst):
        out = []
        if rowfirst:  # MO coefficients: loops imo (outer) iao (inner)
            for imo in range(M.shape[1]):
                for iao in range(M.shape[0]):
                    v = M[iao, imo]
                    if abs(v) >= thresh:
                        out.append(f"{imo + 1:6d}{iao + 1:6d} {v:{fw}.{ndigits}E}")
        else:
            for i in range(M.shape[0]):
                for j in range(M.shape[1]):
                    v = M[i, j]
                    if abs(v) >= thresh:
                        out.append(f"{i + 1:6d}{j + 1:6d} {v:{fw}.{ndigits}E}")
        return out

    for ik in range(nk):
        for sp in range(ns):
            lines.append(f"# BEGIN_KPOINT_SPIN ikp ispin ={ik + 1:6d}{sp + 1:6d}")
            lines.append("# EIGENVALUES")
            lines += vec(eig[(ik, sp)])
            lines.append("# OCCUPATIONS")
            lines += vec(occ[(ik, sp)])
            C = coef[(ik, sp)]
            lines.append("# MO_COEFF_RE (Sparse: mo_index  ao_index  value)")
            lines += sparse(C.real, True)
            if not real:
                lines.append("# MO_COEFF_IM (Sparse: mo_index  ao_index  value)")
                lines += sparse(C.imag, True)
            lines.append(f"# END_KPOINT_SPIN ikp ispin ={ik + 1:6d}{sp + 1:6d}")
        _, S = model.hk(kpts[ik], 0)
        lines.append(f"# BEGIN_OVERLAP ikp ={ik + 1:6d}")
        lines.append("# OVERLAP_RE (Sparse: ao_index_1  ao_index_2  value)")
        lines += sparse(S.real, False)
        if not real:
            lines.append("# OVERLAP_IM (Sparse: ao_index_1  ao_index_2  value)")
            lines += sparse(S.imag, False)
        lines.append(f"# END_OVERLAP ikp ={ik + 1:6d}")
    lines.append("# END_OF_FILE")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_csr(folder, model: Model, kpts, prefix="cohp_ks", suffix="-1_0", real=False):
    folder = Path(folder)
    ns = 2 if model.uks else 1
    for ik, k in enumerate(kpts):
        for sp in range(ns):
            H, _ = model.hk(k, sp)
            lines = []
            for i in range(model.nao):
                for j in range(model.nao):
                    v = H[i, j]
                    if real:
                        lines.append(f"{i + 1:8d}{j + 1:8d}{v.real:23.14E}")
                    else:
                        lines.append(f"{i + 1:8d}{j + 1:8d}{v.real:23.14E}{v.imag:23.14E}")
            (folder / f"{prefix}-KS_SPIN_{sp + 1}_K_{ik + 1}{suffix}.csr").write_text("\n".join(lines) + "\n")


def generate_demo(folder, uks=False, metal=False, mesh=(3, 3, 3), gamma_centered=True, seed=20261122):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    model = Model(uks=uks, metal=metal, seed=seed)
    kpts = mp_mesh(mesh, gamma_centered)
    eig, coef, occ, ef = solve(model, kpts, kT=0.01 if metal else 0.0)
    write_mokp(folder / "cohp_mo-1_0.mokp", model, kpts, eig, coef, occ)
    write_csr(folder, model, kpts)
    (folder / "README_DEMO.txt").write_text(
        "合成紧束缚模型（非 CP2K 计算结果），仅用于试用软件和自检。\n"
        f"原子: {' '.join(model.elements)}；k 网格 {mesh}；{'UKS' if uks else 'RKS'}；{'金属/展宽占据' if metal else '绝缘体'}\n"
        f"模型费米能 {ef * HARTREE_EV:.6f} eV (绝对)\n", encoding="utf-8")
    return folder, model, kpts, ef
