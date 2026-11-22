"""Independent Γ-point supercell reference for the synthetic model."""
import itertools
import numpy as np
from scipy.linalg import eigh


def supercell_matrices(model, N, sp=0):
    nao = model.nao
    cells = list(itertools.product(range(N[0]), range(N[1]), range(N[2])))
    idx = {c: i for i, c in enumerate(cells)}
    n = len(cells) * nao
    H = np.zeros((n, n)); S = np.zeros((n, n))
    for c in cells:
        for R in model.S_R:
            c2 = tuple((c[i] + R[i]) % N[i] for i in range(3))
            i0, j0 = idx[c] * nao, idx[c2] * nao
            S[i0:i0 + nao, j0:j0 + nao] += model.S_R[R]
            H[i0:i0 + nao, j0:j0 + nao] += model.H_R[sp][R]
    return H, S, cells, idx


def reference_integrals(model, N, a, b, r_wrap):
    """ICOHP (eV), ICOOP, ICOBI(Löwdin/Wiberg), ICOBI(Mayer) for bond a(0)-b(R) summed over spins."""
    HARTREE = 27.211386245988
    ns = 2 if model.uks else 1
    g = 2.0 if ns == 1 else 1.0
    sl = model.atom_slices()
    tot = dict(cohp=0.0, coop=0.0, cobi_l=0.0, cobi_m=0.0, cohp_l=0.0)
    ncell = N[0] * N[1] * N[2]
    nel = model.nelectrons() * ncell
    for sp in range(ns):
        H, S, cells, idx = supercell_matrices(model, N, sp)
        e, C = eigh(H, S)
        nfill = int(round(nel / g)) if ns == 1 else None
        # pooled filling across spins for UKS handled by caller (insulator only, same as k-space)
        f = np.zeros(e.size)
        if ns == 1:
            f[:nfill] = g
        else:
            raise NotImplementedError
        c2 = tuple(r % N[i] for i, r in enumerate(r_wrap))
        ia = np.arange(*sl[a - 1]) + idx[(0, 0, 0)] * model.nao
        ib = np.arange(*sl[b - 1]) + idx[c2] * model.nao
        P = (C * f) @ C.T
        tot["cohp"] += np.sum(P[np.ix_(ib, ia)].T * H[np.ix_(ia, ib)]) * 2 * HARTREE
        tot["coop"] += np.sum(P[np.ix_(ib, ia)].T * S[np.ix_(ia, ib)]) * 2
        w, U = np.linalg.eigh(S)
        sh = (U * np.sqrt(w)) @ U.T
        sih = (U / np.sqrt(w)) @ U.T
        PL = sh @ P @ sh
        HL = sih @ H @ sih
        tot["cohp_l"] += np.sum(PL[np.ix_(ib, ia)].T * HL[np.ix_(ia, ib)]) * 2 * HARTREE
        tot["cobi_l"] += np.sum(PL[np.ix_(ia, ib)] ** 2)
        PS = P @ S
        tot["cobi_m"] += np.sum(PS[np.ix_(ia, ib)] * PS[np.ix_(ib, ia)].T)
    return tot
