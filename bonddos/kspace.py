"""k-point system: mesh inference, CP2K gauge, structure helpers, per-k matrices and validation."""
from __future__ import annotations

import itertools
import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .constants import HARTREE_EV
from .mokp import MOKPData, read_mokp, find_ks_sets, read_csr, MOKPError


# ------------------------------------------------------------------ mesh

@dataclass
class MeshInfo:
    mesh: tuple
    shift: tuple
    full: bool
    equal_weights: bool
    messages: list = field(default_factory=list)


def infer_mesh(kpoints: np.ndarray, weights: np.ndarray, tol=1e-5) -> MeshInfo:
    """Infer a regular (possibly shifted) Monkhorst-Pack mesh from a k-point list."""
    msgs = []
    k = np.mod(np.asarray(kpoints, float) + 1e-9, 1.0)
    mesh, shift = [], []
    for ax in range(3):
        vals = np.sort(k[:, ax])
        uniq = [vals[0]]
        for v in vals[1:]:
            if v - uniq[-1] > tol:
                uniq.append(v)
        if len(uniq) > 1 and (1.0 - uniq[-1] + uniq[0]) < tol:
            uniq = uniq[:-1]
        n = len(uniq)
        mesh.append(n)
        sh = uniq[0] * n
        shift.append(round(sh, 6))
        if n > 1:
            steps = np.diff(np.array(uniq + [uniq[0] + 1.0]))
            if np.max(np.abs(steps - 1.0 / n)) > 10 * tol:
                msgs.append(f"第 {ax + 1} 方向 k 点不是等间距网格")
    nk = len(k)
    full = (mesh[0] * mesh[1] * mesh[2] == nk)
    if not full:
        msgs.append(f"k 点数 {nk} ≠ 网格 {mesh[0]}×{mesh[1]}×{mesh[2]}：可能使用了对称性约化（需要 FULL_GRID T, SYMMETRY F）")
    else:
        keys = set(tuple(np.round(row * np.array(mesh) - np.array(shift)).astype(int) % np.array(mesh)) for row in k)
        if len(keys) != nk:
            full = False
            msgs.append("k 点有重复或不在规则网格上")
    eq = bool(np.allclose(weights, weights.mean(), rtol=1e-6, atol=1e-10))
    if not eq:
        msgs.append("k 点权重不相等：网格经过了对称性约化，周期镜像（IMAGE）分析不可用")
    if abs(weights.sum() - 1.0) > 1e-6:
        msgs.append(f"k 点权重之和为 {weights.sum():.8f}（应为 1）")
    return MeshInfo(tuple(mesh), tuple(shift), full and eq, eq, msgs)


# ------------------------------------------------------------------ structure helpers

def anint(x):
    """Fortran ANINT: round half away from zero."""
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def frac_coords(cell, pos):
    return np.linalg.solve(cell.T, np.asarray(pos).T).T


def parse_periodic(text: str):
    text = (text or "XYZ").upper().replace(" ", "")
    if text in ("NONE", "0"):
        return (False, False, False)
    return tuple(ax in text for ax in "XYZ")


@dataclass
class Bond:
    a: int                  # 1-based atom index (home cell)
    b: int                  # 1-based atom index
    r_raw: tuple            # lattice translation in the raw (printed) coordinate frame
    r_wrap: tuple           # translation in CP2K's internal wrapped frame (used for H(R))
    length: float           # Angstrom
    mode: str = "IMAGE"     # IMAGE | PAIR
    label: str = ""

    def key(self):
        return (self.a, self.b, tuple(self.r_raw), self.mode)


class Structure:
    def __init__(self, cell, atoms, periodic=(True, True, True)):
        self.cell = np.asarray(cell, float)
        self.atoms = atoms
        self.periodic = tuple(periodic)
        self.frac_raw = np.array([frac_coords(self.cell, a.pos) for a in atoms])
        per = np.array(self.periodic, bool)
        g = np.zeros_like(self.frac_raw)
        g[:, per] = -anint(self.frac_raw[:, per])
        self.g = g.astype(int)                     # s_wrapped = s_raw + g
        self.frac_wrap = self.frac_raw + self.g
        inv_norm = np.linalg.norm(np.linalg.inv(self.cell), axis=0)   # |b_i|/2pi per axis
        tol = 2.0e-6 * inv_norm * math.sqrt(3)
        frac_part = np.abs(np.abs(self.frac_raw) - np.floor(np.abs(self.frac_raw)) - 0.5)
        self.ambiguous = [(i + 1, ax) for i in range(len(atoms)) for ax in range(3)
                          if per[ax] and frac_part[i, ax] < tol[ax]]

    def symbol(self, i):
        return self.atoms[i - 1].element

    def label(self, i):
        return f"{self.atoms[i - 1].element}{i}"

    def wrap_translation(self, a, b, r_raw):
        """R_wrap = R_raw - g_B + g_A."""
        return tuple(int(x) for x in np.asarray(r_raw) - self.g[b - 1] + self.g[a - 1])

    def raw_translation(self, a, b, r_wrap):
        return tuple(int(x) for x in np.asarray(r_wrap) + self.g[b - 1] - self.g[a - 1])

    def bond_vector(self, a, b, r_raw):
        ds = self.frac_raw[b - 1] + np.asarray(r_raw, float) - self.frac_raw[a - 1]
        return ds @ self.cell

    def distance(self, a, b, r_raw):
        return float(np.linalg.norm(self.bond_vector(a, b, r_raw)))

    def _search_range(self, dmax):
        inv = np.linalg.inv(self.cell)          # columns = b_i / 2pi
        n = [int(math.ceil(dmax * np.linalg.norm(inv[:, i]))) + 1 if self.periodic[i] else 0 for i in range(3)]
        return n

    def images(self, a, b, dmax):
        """All (R_raw, d) with d <= dmax for the pair (a, b)."""
        n = self._search_range(dmax)
        ds0 = self.frac_raw[b - 1] - self.frac_raw[a - 1]
        base = [-int(round(ds0[i])) if self.periodic[i] else 0 for i in range(3)]
        rng = [range(base[i] - n[i], base[i] + n[i] + 1) for i in range(3)]
        out = []
        for R in itertools.product(*rng):
            v = (ds0 + np.array(R)) @ self.cell
            d = float(np.linalg.norm(v))
            if d <= dmax + 1e-9:
                out.append((R, d))
        out.sort(key=lambda t: (round(t[1], 6), [abs(x) for x in t[0]], t[0]))
        return out

    def nearest(self, a, b, all_equivalent=False, tol=1e-5):
        dmax = 1.0
        while True:
            im = [(R, d) for R, d in self.images(a, b, dmax) if not (a == b and R == (0, 0, 0))]
            if im:
                break
            dmax *= 2
            if dmax > 1e4:
                raise ValueError("找不到最近镜像")
        dmin = im[0][1]
        # make sure no shorter image exists beyond the search box
        im = [(R, d) for R, d in self.images(a, b, dmin + 1e-3) if not (a == b and R == (0, 0, 0))]
        dmin = im[0][1]
        sel = [(R, d) for R, d in im if d <= dmin + tol]
        if a == b:
            sel = _dedupe_self(sel)
        return sel if all_equivalent else sel[:1]

    def neighbor_search(self, dmin, dmax, elements_a=None, elements_b=None, atoms_subset=None):
        """All symmetry-unique bonds (a<=b, dedup R/-R for a==b) with dmin <= d <= dmax."""
        out = []
        n = len(self.atoms)
        subset = set(atoms_subset) if atoms_subset else None
        for a in range(1, n + 1):
            for b in range(a, n + 1):
                if subset and a not in subset and b not in subset:
                    continue
                ea, eb = self.symbol(a), self.symbol(b)
                if elements_a and elements_b:
                    ok = (ea in elements_a and eb in elements_b) or (eb in elements_a and ea in elements_b)
                    if not ok:
                        continue
                elif elements_a:
                    if ea not in elements_a and eb not in elements_a:
                        continue
                ims = self.images(a, b, dmax)
                if a == b:
                    ims = _dedupe_self([(R, d) for R, d in ims if R != (0, 0, 0)])
                for R, d in ims:
                    if d >= dmin - 1e-9:
                        out.append((a, b, R, d))
        out.sort(key=lambda t: (round(t[3], 5), t[0], t[1], t[2]))
        return out


def _dedupe_self(items):
    seen = set()
    out = []
    for R, d in items:
        neg = tuple(-x for x in R)
        if neg in seen:
            continue
        seen.add(tuple(R))
        out.append((R, d))
    return out


# ------------------------------------------------------------------ system

@dataclass
class CheckResult:
    name: str
    value: float
    ok: bool
    note: str = ""


class KSystem:
    """CP2K k-point data (C, S, H, eigenvalues, occupations) with validation and Löwdin helpers."""

    def __init__(self, mokp: MOKPData, ks_files: dict | None = None, triangle="FULL", periodic="XYZ",
                 progress=None):
        self.d = mokp
        self.nk = mokp.nkp
        self.ns = mokp.nspins
        self.nao = mokp.nao
        self.nmo = mokp.nmo
        self.g = mokp.degeneracy()
        self.kpts = mokp.kpoints
        self.wts = mokp.weights
        self.mesh = infer_mesh(self.kpts, self.wts)
        self.structure = Structure(mokp.cell, mokp.atoms, parse_periodic(periodic))
        self.H = {}
        self.h_source = "none"
        self.warnings = list(mokp.warnings)
        self.checks: list[CheckResult] = []
        self._lowdin = {}
        if ks_files:
            need = [(s, k) for s in range(self.ns) for k in range(self.nk)]
            missing = [(s + 1, k + 1) for s, k in need if (s, k) not in ks_files]
            if missing:
                raise MOKPError(f"KS 文件缺少 (spin,k): {missing[:6]}{' ...' if len(missing) > 6 else ''}")
            for i, (s, k) in enumerate(need):
                if progress:
                    progress(i / len(need), f"读取 KS 矩阵 {i + 1}/{len(need)}")
                self.H[(k, s)] = read_csr(ks_files[(s, k)], self.nao, triangle)
            self.h_source = "KS_CSR_WRITE"
        elif self.nmo == self.nao and mokp.has_overlap:
            for k in range(self.nk):
                S = mokp.overlap[k]
                for s in range(self.ns):
                    C = mokp.coeff[(k, s)]
                    SC = S @ C
                    self.H[(k, s)] = (SC * mokp.eig[(k, s)][None, :]) @ SC.conj().T
            self.h_source = "reconstructed (S C ε C† S)"
            self.warnings.append("未提供 KS 文件，H(k) 由 S·C·ε·C†·S 重建（仅当 NMO=NAO 时严格成立）")
        else:
            self.warnings.append("未提供 KS_CSR 文件且 NMO<NAO：COHP 不可用（COOP/COBI/电荷/PDOS 仍可计算）")
        if self.structure.ambiguous:
            self.warnings.append("有原子分数坐标恰好位于 ±0.5 边界附近，CP2K 内部回卷方向可能不确定："
                                 + ", ".join(f"{self.structure.label(i)}({'abc'[ax]})" for i, ax in self.structure.ambiguous[:8])
                                 + "。周期镜像标签已自动检验（见报告），如仍有疑问可将整体结构平移 0.01 Å 重算")

    # ---- basic accessors
    @property
    def has_H(self):
        return bool(self.H)

    @property
    def has_S(self):
        return self.d.has_overlap

    def S(self, k):
        return self.d.overlap[k]

    def C(self, k, s):
        return self.d.coeff[(k, s)]

    def eps_ev(self, k, s):
        return self.d.eig[(k, s)] * HARTREE_EV

    def occ(self, k, s):
        return self.d.occ[(k, s)]

    def H_ev(self, k, s):
        return self.H[(k, s)] * HARTREE_EV

    def lowdin(self, k):
        """(S^{1/2}, S^{-1/2}) for k (cached)."""
        if k not in self._lowdin:
            S = self.S(k)
            S = 0.5 * (S + S.conj().T)
            w, U = np.linalg.eigh(S)
            if w.min() <= 0:
                raise MOKPError(f"k 点 {k + 1} 的重叠矩阵不是正定的 (最小本征值 {w.min():.3e})")
            sh = (U * np.sqrt(w)) @ U.conj().T
            sih = (U / np.sqrt(w)) @ U.conj().T
            self._lowdin[k] = (sh, sih, w.min(), w.max())
        return self._lowdin[k]

    def atom_aos(self, i):
        a = self.d.atoms[i - 1]
        return np.arange(a.first_ao - 1, a.last_ao)

    # ---- reference energies
    def band_edges(self):
        occ_e, vir_e = [], []
        partial = False
        for (k, s), o in self.d.occ.items():
            e = self.eps_ev(k, s)
            rel = o / self.g
            occ_e.append(e[rel > 0.5])
            vir_e.append(e[rel <= 0.5])
            if np.any((rel > 1e-6) & (rel < 1 - 1e-6)):
                partial = True
        occ_e = np.concatenate(occ_e) if occ_e else np.array([])
        vir_e = np.concatenate(vir_e) if vir_e else np.array([])
        vbm = float(occ_e.max()) if occ_e.size else None
        cbm = float(vir_e.min()) if vir_e.size else None
        return vbm, cbm, partial

    def fermi_estimate(self):
        """Chemical potential estimate: energy where the occupation fraction crosses 1/2."""
        es, fs = [], []
        for (k, s), o in self.d.occ.items():
            es.append(self.eps_ev(k, s))
            fs.append(o / self.g)
        e = np.concatenate(es)
        f = np.concatenate(fs)
        order = np.argsort(e)
        e, f = e[order], f[order]
        idx = np.where(f < 0.5)[0]
        if idx.size == 0:
            return float(e[-1])
        i = idx[0]
        if i == 0:
            return float(e[0])
        f0, f1 = f[i - 1], f[i]
        if f0 == f1:
            return float(0.5 * (e[i - 1] + e[i]))
        t = (f0 - 0.5) / (f0 - f1)
        return float(e[i - 1] + t * (e[i] - e[i - 1]))

    def reference_energy(self, mode="auto", value=None):
        """Absolute eV used as energy zero. mode: auto|hoco|midgap|fermi|manual."""
        vbm, cbm, partial = self.band_edges()
        mode = (mode or "auto").lower()
        if mode == "manual":
            return float(value), "手动指定"
        if mode == "fermi" or (mode == "auto" and partial):
            return self.fermi_estimate(), "由占据数估计的费米能 (f=1/2)"
        if mode == "midgap" and vbm is not None and cbm is not None:
            return 0.5 * (vbm + cbm), "带隙中点"
        if vbm is not None:
            return vbm, "最高占据态 (HOCO/VBM)"
        return self.fermi_estimate(), "费米能估计"

    # ---- validation
    def validate(self, progress=None, electrons=None):
        checks = []
        max_herm = max_norm = max_res = 0.0
        min_seig = np.inf
        max_cond = 0.0
        trace_diff = 0.0
        band_energy = 0.0
        for k in range(self.nk):
            if progress:
                progress(k / self.nk, f"数据检查 k {k + 1}/{self.nk}")
            if self.has_S:
                S = self.S(k)
                max_herm = max(max_herm, float(np.max(np.abs(S - S.conj().T))))
                w = np.linalg.eigvalsh(0.5 * (S + S.conj().T))
                min_seig = min(min_seig, float(w.min()))
                max_cond = max(max_cond, float(w.max() / max(w.min(), 1e-300)))
            for s in range(self.ns):
                C = self.C(k, s)
                e = self.d.eig[(k, s)]
                f = self.occ(k, s)
                band_energy += self.wts[k] * float(np.dot(f, e)) * HARTREE_EV
                if self.has_S:
                    SC = self.S(k) @ C
                    ov = C.conj().T @ SC
                    max_norm = max(max_norm, float(np.max(np.abs(ov - np.eye(self.nmo)))))
                    if self.has_H and self.h_source == "KS_CSR_WRITE":
                        Hk = self.H[(k, s)]
                        r = Hk @ C - SC * e[None, :]
                        scale = max(1.0, float(np.max(np.abs(e))))
                        max_res = max(max_res, float(np.max(np.abs(r))) / scale)
                        Pk = (C * f[None, :]) @ C.conj().T
                        trace_diff += self.wts[k] * float(np.real(np.trace(Pk @ Hk)) - np.dot(f, e)) * HARTREE_EV
        if self.has_S:
            checks.append(CheckResult("S(k) 厄米性 max|S−S†|", max_herm, max_herm < 1e-6))
            checks.append(CheckResult("S(k) 最小本征值", min_seig, min_seig > 1e-10,
                                      "接近 0 表示基组近线性相关，Löwdin 结果需谨慎"))
            checks.append(CheckResult("S(k) 条件数", max_cond, max_cond < 1e10))
            checks.append(CheckResult("MO 归一化 max|C†SC−I|", max_norm, max_norm < 1e-4,
                                      "偏大时检查 MO_KP 的 NDIGITS（建议 14）"))
        if self.h_source == "KS_CSR_WRITE":
            checks.append(CheckResult("广义本征方程 max|HC−SCε|/max|ε| (Ha)", max_res, max_res < 1e-4,
                                      "偏大说明 KS 文件与 MO_KP 不是同一电子步/同一 k 点顺序"))
            checks.append(CheckResult("Tr(PH) − Σwfε (eV)", trace_diff, abs(trace_diff) < 1e-3))
        nel = self.d.electrons()
        ok = True if electrons is None else abs(nel - electrons) < 1e-4
        checks.append(CheckResult("电子数 Σ w·f", nel, ok, "" if electrons is None else f"期望 {electrons}"))
        # last exported MO must be (almost) empty
        last_occ = max(float(self.occ(k, s)[-1]) for k in range(self.nk) for s in range(self.ns))
        checks.append(CheckResult("最高导出轨道占据数", last_occ, last_occ < 1e-6,
                                  "过大说明 ADDED_MOS 太少，可能遗漏占据态"))
        if self.has_S and self.mesh.full and self.nk > 1:
            # S(R=0) onsite diagonal should be exactly 1 for normalized AOs; deviation = aliasing of images
            dev = self.onsite_alias_error()
            checks.append(CheckResult("k 网格混叠指标 max|S̃_μμ(0)−1|", dev, dev < 1e-2,
                                      "反映有限 k 网格把远处镜像折叠回来的程度；IMAGE 模式需要它足够小"))
        self.checks = checks
        self.band_energy = band_energy
        if progress:
            progress(1.0, "检查完成")
        return checks

    def onsite_alias_error(self):
        acc = np.zeros(self.nao, dtype=complex)
        for k in range(self.nk):
            acc += self.wts[k] * np.diag(self.S(k))
        return float(np.max(np.abs(acc - 1.0)))

    # ---- image / aliasing helpers
    def image_status(self, r_wrap):
        """Classify a wrapped translation against the k mesh: ok | nyquist | aliased | gamma."""
        status = "ok"
        notes = []
        for ax in range(3):
            n = self.mesh.mesh[ax]
            r = r_wrap[ax]
            if not self.structure.periodic[ax]:
                if r != 0:
                    return "nonperiodic", [f"{'abc'[ax]} 方向非周期"]
                continue
            if n == 1:
                status = "gamma"
                notes.append(f"{'abc'[ax]} 方向只有 1 个 k 点：该方向所有周期镜像叠加在一起（表面/二维体系请把周期方向设为 XY 等）")
                continue
            if 2 * abs(r) > n:
                return "aliased", [f"{'abc'[ax]} 方向 |R|={abs(r)} 超过网格 {n} 的一半"]
            if 2 * abs(r) == n:
                status = "nyquist" if status == "ok" else status
                notes.append(f"{'abc'[ax]} 方向 R={r} 恰在 Nyquist 边界（R 与 −R 混叠）")
        return status, notes


def load_system(mokp_path, ks_dir=None, ks_set=None, triangle="FULL", periodic="XYZ", progress=None):
    """Convenience loader. ks_set: label from find_ks_sets (None -> the only set / none)."""
    def p_mokp(f, t):
        if progress:
            progress(0.6 * f, t)
    data = read_mokp(mokp_path, progress=p_mokp)
    files = None
    if ks_dir:
        sets = find_ks_sets(ks_dir)
        if sets:
            if ks_set is None:
                if len(sets) > 1:
                    # prefer the set that covers all (spin,k)
                    full = [lab for lab, fs in sets.items() if len(fs) == data.nkp * data.nspins]
                    if len(full) == 1:
                        ks_set = full[0]
                    else:
                        raise MOKPError("KS 文件夹中有多组文件，请选择一组: " + "; ".join(sets))
                else:
                    ks_set = next(iter(sets))
            files = sets[ks_set]

    def p_ks(f, t):
        if progress:
            progress(0.6 + 0.4 * f, t)
    return KSystem(data, files, triangle=triangle, periodic=periodic, progress=p_ks)


def parse_atom_selection(text: str, natom: int, structure: Structure | None = None):
    """'1-4,7' or element symbols 'W,O' or 'all' -> sorted list of 1-based atom indices."""
    text = (text or "").strip()
    if not text or text.lower() in ("all", "全部", "*"):
        return list(range(1, natom + 1))
    out = set()
    for tok in re.split(r"[,\s;，]+", text):
        if not tok:
            continue
        m = re.fullmatch(r"(\d+)\s*-\s*(\d+)", tok)
        if m:
            a, b = int(m.group(1)), int(m.group(2))
            out.update(range(min(a, b), max(a, b) + 1))
        elif tok.isdigit():
            out.add(int(tok))
        elif structure is not None:
            hits = [i for i in range(1, natom + 1) if structure.symbol(i).lower() == tok.lower()]
            if not hits:
                raise ValueError(f"未知的原子/元素: {tok}")
            out.update(hits)
        else:
            raise ValueError(f"无法解析原子选择: {tok}")
    bad = [i for i in out if i < 1 or i > natom]
    if bad:
        raise ValueError(f"原子编号超出范围 1..{natom}: {bad[:5]}")
    return sorted(out)
