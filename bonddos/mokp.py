"""Reader for CP2K 2026.2 ``&MO_KP`` (``.mokp``, CP2K_KPOINT_MO_DUMP Version 2.0) and KS CSR files.

Format facts (CP2K v2026.2 ``kpoint_mo_dump.F`` / DBCSR ``csr_write``)
* energies always Hartree; lengths Angstrom or Bohr (header); k in fractional reciprocal coordinates
* ATOM_LIST rows: id, element, Z_nuc, Z_eff, x, y, z, first_AO, last_AO
* AO_LIST rows (OVERLAP_MATRIX mode): ao, atom, atom_shell_index, l, m   (m = -l..l)
* sparse MO records: ``mo_index ao_index value`` with format (2I6,1X,ES..)
* S(k) is written once per k (independent of spin); imaginary parts only for complex wfn
* KS CSR (k-space, REAL_SPACE F): ``row col re [im]`` (2I8, 2E23.14), 1-based, no header
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .constants import BOHR_ANG


class MOKPError(ValueError):
    pass


@dataclass(eq=False)
class Atom:
    index: int          # 1-based
    element: str
    z_nuc: int
    z_eff: float
    pos: np.ndarray     # Angstrom (raw printed coordinates)
    first_ao: int       # 1-based inclusive
    last_ao: int

    @property
    def label(self):
        return f"{self.element}{self.index}"

    @property
    def ao_slice(self):
        return slice(self.first_ao - 1, self.last_ao)


@dataclass(eq=False)
class MOKPData:
    path: str
    natom: int
    nspins: int
    nao: int
    nkp: int
    nmo: int
    use_real_wfn: bool
    ao_export: str
    cell: np.ndarray                    # 3x3 Angstrom, rows = lattice vectors a, b, c
    atoms: list
    kpoints: np.ndarray                 # (nkp, 3) fractional
    weights: np.ndarray                 # (nkp,)
    ao_atom: np.ndarray                 # (nao,) 1-based atom index
    ao_shell: np.ndarray                # (nao,) shell index within atom (1-based)
    ao_l: np.ndarray
    ao_m: np.ndarray
    eig: dict = field(default_factory=dict)     # (ik, is) -> (nmo,) Hartree, 0-based keys
    occ: dict = field(default_factory=dict)     # (ik, is) -> (nmo,)
    coeff: dict = field(default_factory=dict)   # (ik, is) -> (nao, nmo) complex
    overlap: dict = field(default_factory=dict) # ik -> (nao, nao) complex
    sparse_threshold: float = 0.0
    warnings: list = field(default_factory=list)

    @property
    def has_overlap(self):
        return len(self.overlap) == self.nkp

    def atom_of(self, i):
        return self.atoms[i - 1]

    def electrons(self):
        return float(sum(self.weights[k] * self.occ[(k, s)].sum() for (k, s) in self.occ))

    def degeneracy(self):
        """Spin degeneracy of one channel: 2 for restricted, 1 for UKS."""
        return 2.0 if self.nspins == 1 else 1.0


# Fortran ES/E edit descriptors drop the exponent letter when the exponent needs three digits
# (e.g. Fermi-smeared occupations such as ``2.91755220991368-102``).
_BARE_EXP = re.compile(r"(?<=[0-9.])([+-]\d{2,4})(?![0-9.])")


def _to_floats(text):
    """Whitespace separated Fortran reals -> 1-D float array (handles D exponents and ``1.0-102``)."""
    toks = text.split()
    try:
        return np.array(toks, dtype=float)
    except ValueError:
        text = text.replace("D", "E").replace("d", "e")
        try:
            return np.array(text.split(), dtype=float)
        except ValueError:
            return np.array(_BARE_EXP.sub(r"E\1", text).split(), dtype=float)


def _floats(lines):
    return _to_floats(" ".join(lines))


def _sparse(lines, ncol_expected=3):
    """Parse (2I6,1X,ESw.d) sparse records robustly (indices may fuse above 99999)."""
    if not lines:
        return np.zeros((0, 3))
    toks = " ".join(lines).split()
    if len(toks) == 3 * len(lines):
        return _to_floats(" ".join(toks)).reshape(-1, 3)
    out = np.empty((len(lines), 3))
    for i, ln in enumerate(lines):
        p = ln.split()
        if len(p) == 3:
            out[i] = [float(p[0]), float(p[1]), _to_floats(p[2])[0]]
        else:  # fixed width fallback
            out[i] = [float(ln[0:6]), float(ln[6:12]), _to_floats(ln[12:])[0]]
    return out


def read_mokp(path, progress=None, keep_coefficients=True) -> MOKPData:
    """Parse a .mokp file. ``progress(fraction, text)`` is called occasionally."""
    path = Path(path)
    size = max(1, path.stat().st_size)
    header = {}
    cell_rows = []
    atoms = []
    kpts = []
    ao_rows = []
    blocks = {}
    overlaps = {}
    section = None
    buf = []
    cur_k = cur_s = None
    sub = None
    seen_end = False
    length_unit = "bohr"
    read_bytes = 0
    last_report = 0.0
    nao = None
    nmo = None
    tmp = {}

    def flush_sub():
        nonlocal buf, sub
        if sub is None:
            buf = []
            return
        if section == "kspin":
            if sub == "eig":
                tmp["eig"] = _floats(buf)
            elif sub == "occ":
                tmp["occ"] = _floats(buf)
            elif sub == "cre":
                tmp["cre"] = _sparse(buf)
            elif sub == "cim":
                tmp["cim"] = _sparse(buf)
        elif section == "ovl":
            if sub == "sre":
                tmp["sre"] = _sparse(buf)
            elif sub == "sim":
                tmp["sim"] = _sparse(buf)
        buf = []
        sub = None

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            read_bytes += len(line)
            if progress and read_bytes / size - last_report > 0.02:
                last_report = read_bytes / size
                progress(min(0.99, last_report), "读取 MO_KP 文件")
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                t = s[1:].strip()
                tu = t.upper()
                if tu.startswith("CP2K_KPOINT_MO_DUMP"):
                    header["version"] = t
                    continue
                if tu.startswith("ALL ENERGIES"):
                    length_unit = "angstrom" if "ANGSTROM" in tu else "bohr"
                    continue
                if tu.startswith("DIMENSIONS"):
                    nums = [int(x) for x in re.findall(r"-?\d+", t.split("=", 1)[1])]
                    header["natom"], header["nspins"], header["nao"], header["nkp"] = nums[:4]
                    nao = header["nao"]
                    continue
                if tu.startswith("NMO"):
                    header["nmo"] = int(t.split("=", 1)[1])
                    nmo = header["nmo"]
                    continue
                if tu.startswith("USE_REAL_WFN"):
                    header["use_real"] = "T" in t.split("=", 1)[1].upper()
                    continue
                if tu.startswith("AO_EXPORT_TYPE"):
                    header["ao_export"] = t.split("=", 1)[1].strip().upper()
                    continue
                if tu.startswith("SPARSE_THRESHOLD"):
                    header["thresh"] = float(t.split("=", 1)[1].replace("D", "E"))
                    continue
                if tu.startswith("CELL_VECTORS"):
                    section = "cell"
                    continue
                if tu.startswith("ATOM_LIST"):
                    section = "atoms"
                    continue
                if tu.startswith("KPOINT_LIST"):
                    section = "kpts"
                    continue
                if tu.startswith("AO_LIST"):
                    section = "aolist"
                    continue
                if tu.startswith("GTO_BASIS"):
                    section = "gto"
                    continue
                if tu.startswith("BEGIN_KPOINT_SPIN"):
                    nums = [int(x) for x in re.findall(r"\d+", t.split("=", 1)[1])]
                    cur_k, cur_s = nums[0] - 1, nums[1] - 1
                    section = "kspin"
                    tmp.clear()
                    continue
                if tu.startswith("END_KPOINT_SPIN"):
                    flush_sub()
                    if "eig" not in tmp or "occ" not in tmp:
                        raise MOKPError(f"k 点 {cur_k + 1} 自旋 {cur_s + 1} 缺少本征值或占据数")
                    key = (cur_k, cur_s)
                    if key in blocks:
                        raise MOKPError(f"k 点 {cur_k + 1} 自旋 {cur_s + 1} 重复出现：文件可能包含多次输出，请只保留最后一次")
                    if keep_coefficients:
                        c = np.zeros((nao, nmo), dtype=complex)
                        cre = tmp.get("cre", np.zeros((0, 3)))
                        if cre.size:
                            mo = cre[:, 0].astype(int) - 1
                            ao = cre[:, 1].astype(int) - 1
                            c[ao, mo] = cre[:, 2]
                        cim = tmp.get("cim")
                        if cim is not None and cim.size:
                            mo = cim[:, 0].astype(int) - 1
                            ao = cim[:, 1].astype(int) - 1
                            c[ao, mo] += 1j * cim[:, 2]
                    else:
                        c = None
                    blocks[key] = (tmp["eig"][:nmo], tmp["occ"][:nmo], c)
                    tmp.clear()
                    section = None
                    continue
                if tu.startswith("BEGIN_OVERLAP"):
                    cur_k = int(re.findall(r"\d+", t.split("=", 1)[1])[0]) - 1
                    section = "ovl"
                    tmp.clear()
                    continue
                if tu.startswith("END_OVERLAP"):
                    flush_sub()
                    smat = np.zeros((nao, nao), dtype=complex)
                    sre = tmp.get("sre", np.zeros((0, 3)))
                    if sre.size:
                        smat[sre[:, 0].astype(int) - 1, sre[:, 1].astype(int) - 1] = sre[:, 2]
                    sim = tmp.get("sim")
                    if sim is not None and sim.size:
                        smat[sim[:, 0].astype(int) - 1, sim[:, 1].astype(int) - 1] += 1j * sim[:, 2]
                    overlaps[cur_k] = smat
                    tmp.clear()
                    section = None
                    continue
                if tu.startswith("EIGENVALUES"):
                    flush_sub(); sub = "eig"; continue
                if tu.startswith("OCCUPATIONS"):
                    flush_sub(); sub = "occ"; continue
                if tu.startswith("MO_COEFF_RE"):
                    flush_sub(); sub = "cre"; continue
                if tu.startswith("MO_COEFF_IM"):
                    flush_sub(); sub = "cim"; continue
                if tu.startswith("OVERLAP_RE"):
                    flush_sub(); sub = "sre"; continue
                if tu.startswith("OVERLAP_IM"):
                    flush_sub(); sub = "sim"; continue
                if tu.startswith("END_OF_FILE"):
                    seen_end = True
                    section = None
                    continue
                continue
            # data line
            if section == "cell":
                cell_rows.append([float(x) for x in s.split()[:3]])
                if len(cell_rows) == 3:
                    section = None
            elif section == "atoms":
                p = s.split()
                atoms.append(p)
            elif section == "kpts":
                p = s.split()
                kpts.append([float(x) for x in p[1:5]])
            elif section == "aolist":
                ao_rows.append([int(x) for x in s.split()[:5]])
            elif section in ("kspin", "ovl") and sub is not None:
                buf.append(s)
    if "nao" not in header:
        raise MOKPError("不是有效的 CP2K MO_KP 文件（缺少 DIMENSIONS 行）")
    if not seen_end:
        raise MOKPError("MO_KP 文件不完整（缺少 END_OF_FILE），计算可能未正常结束或文件被截断")
    scale = 1.0 if length_unit == "angstrom" else BOHR_ANG
    cell = np.array(cell_rows, float) * scale
    at = []
    for p in atoms:
        at.append(Atom(index=int(p[0]), element=p[1], z_nuc=int(p[2]), z_eff=float(p[3]),
                       pos=np.array([float(p[4]), float(p[5]), float(p[6])]) * scale,
                       first_ao=int(p[7]), last_ao=int(p[8])))
    if len(at) != header["natom"]:
        raise MOKPError(f"ATOM_LIST 行数 {len(at)} 与 natom={header['natom']} 不符")
    kp = np.array(kpts, float)
    if kp.shape[0] != header["nkp"]:
        raise MOKPError(f"KPOINT_LIST 行数 {kp.shape[0]} 与 nkp={header['nkp']} 不符")
    nao = header["nao"]
    ao_atom = np.zeros(nao, int)
    ao_shell = np.zeros(nao, int)
    ao_l = np.zeros(nao, int)
    ao_m = np.zeros(nao, int)
    if ao_rows:
        a = np.array(ao_rows, int)
        idx = a[:, 0] - 1
        ao_atom[idx], ao_shell[idx], ao_l[idx], ao_m[idx] = a[:, 1], a[:, 2], a[:, 3], a[:, 4]
    else:
        for atm in at:
            ao_atom[atm.first_ao - 1:atm.last_ao] = atm.index
    wts = kp[:, 3].copy()
    if wts.size and abs(wts.sum() - 1.0) < 1e-6:
        wts /= wts.sum()                 # weights are printed with 9 digits only
    data = MOKPData(path=str(path), natom=header["natom"], nspins=header["nspins"], nao=nao, nkp=header["nkp"],
                    nmo=header.get("nmo", 0), use_real_wfn=header.get("use_real", False),
                    ao_export=header.get("ao_export", "?"), cell=cell, atoms=at, kpoints=kp[:, :3],
                    weights=wts, ao_atom=ao_atom, ao_shell=ao_shell, ao_l=ao_l, ao_m=ao_m,
                    sparse_threshold=header.get("thresh", 0.0))
    for (k, s), (e, o, c) in blocks.items():
        data.eig[(k, s)] = e
        data.occ[(k, s)] = o
        if c is not None:
            data.coeff[(k, s)] = c
    data.overlap = overlaps
    missing = [(k + 1, s + 1) for k in range(data.nkp) for s in range(data.nspins) if (k, s) not in data.eig]
    if missing:
        raise MOKPError(f"缺少 k/自旋数据块: {missing[:6]}{' ...' if len(missing) > 6 else ''}")
    if data.ao_export != "OVERLAP_MATRIX":
        data.warnings.append("MO_KP 未使用 AO_EXPORT_TYPE OVERLAP_MATRIX：没有 S(k)，COOP/COBI/Mulliken 和 PDOS 不可用")
    elif not data.has_overlap:
        raise MOKPError("AO_EXPORT_TYPE=OVERLAP_MATRIX 但缺少部分 k 点的重叠矩阵")
    if progress:
        progress(1.0, "MO_KP 读取完成")
    return data


# ---------------------------------------------------------------- KS CSR files

_KS_RE = re.compile(r"^(?P<prefix>.*?)KS_SPIN_(?P<spin>\d+)_K_(?P<k>\d+)(?P<suffix>.*)\.csr$", re.I)


def find_ks_sets(folder):
    """Group KS CSR files by (prefix, suffix). Returns {set_label: {(ispin0, ik0): Path}}."""
    folder = Path(folder)
    sets = {}
    if not folder.is_dir():
        return sets
    for p in sorted(folder.iterdir()):
        m = _KS_RE.match(p.name)
        if not m:
            continue
        label = f"{m.group('prefix')}KS_SPIN_*_K_*{m.group('suffix')}.csr"
        sets.setdefault(label, {})[(int(m.group("spin")) - 1, int(m.group("k")) - 1)] = p
    return sets


def diagnose_ks_folder(folder):
    """Hints for a folder in which no usable ``*KS_SPIN_<s>_K_<k>*.csr`` file was found."""
    folder = Path(folder)
    if not folder.is_dir():
        return [f"KS 文件夹不存在: {folder}"]
    names = [p.name for p in folder.iterdir() if p.is_file()]
    csr = [n for n in names if n.lower().endswith(".csr")]
    if any(re.search(r"KS_SPIN_\d+-", n, re.I) for n in csr):
        return ["找到的是没有 k 点标记的 KS 文件（…KS_SPIN_1-1_0.csr）：CP2K 输入里没有 &KPOINTS。"
                "只算 Γ 点也必须写 &KPOINTS / SCHEME GAMMA，否则不会输出 .mokp 和逐 k 点的 KS 矩阵"]
    if any(re.search(r"KS_SPIN_\d+_R_", n, re.I) for n in csr):
        return ["找到的是实空间 KS 矩阵（KS_CSR_WRITE 的 REAL_SPACE T）：请改为 REAL_SPACE F 输出 H(k)"]
    if csr:
        return ["文件夹里的 .csr 文件名不含 KS_SPIN_<s>_K_<k>，无法识别为 KS_CSR_WRITE 的 k 空间输出"]
    return ["KS 文件夹里没有 .csr 文件：请在 CP2K 输入的 &DFT/&PRINT 中加入 &KS_CSR_WRITE ON 段"]


def read_csr(path, nao, triangle="FULL") -> np.ndarray:
    """Read a DBCSR text CSR dump into a dense (nao, nao) complex matrix (Hartree)."""
    path = Path(path)
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith(("#", "!"))]
    if not lines:
        raise MOKPError(f"{path.name}: 空文件")
    ncol = len(lines[0].split())
    if ncol not in (3, 4):
        raise MOKPError(f"{path.name}: 每行应为 3 或 4 列（行 列 实部 [虚部]），实际 {ncol} 列；不支持 BINARY T")
    arr = _to_floats(" ".join(lines))
    if arr.size % ncol:
        raise MOKPError(f"{path.name}: 列数不一致，文件可能被截断")
    arr = arr.reshape(-1, ncol)
    r = arr[:, 0].astype(int) - 1
    c = arr[:, 1].astype(int) - 1
    if r.min() < 0 or c.min() < 0 or r.max() >= nao or c.max() >= nao:
        raise MOKPError(f"{path.name}: 矩阵索引超出 AO 数 {nao}，KS 文件与 MO_KP 不是同一体系？")
    m = np.zeros((nao, nao), dtype=complex)
    vals = arr[:, 2] + (1j * arr[:, 3] if ncol == 4 else 0.0)
    m[r, c] = vals
    upper = triangle.upper() == "UPPER"
    if not upper and nao > 1:
        # UPPER_TRIANGULAR T that was not declared: strictly lower triangle empty, upper triangle not
        upper = (not np.any(np.tril(m, -1))) and bool(np.any(np.triu(m, 1)))
    if upper:
        low = np.tril_indices(nao, -1)
        m[low] = np.conj(m.T[low])
    return m
