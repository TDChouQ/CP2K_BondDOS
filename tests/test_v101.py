"""Regression tests added in v1.0.1 (bugs found with real CP2K 2026.2 output)."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from bonddos.constants import HARTREE_EV  # noqa: E402
from bonddos.demo import generate_demo  # noqa: E402
from bonddos.kspace import load_system  # noqa: E402
from bonddos.mokp import _to_floats, read_mokp  # noqa: E402
from bonddos.bonding import BondingEngine, BondingOptions  # noqa: E402
from bonddos.session import BondSession, BondRequest, parse_config, session_from_config  # noqa: E402
from bonddos.spectra import curve_descriptors, make_grid, broaden, trapz  # noqa: E402

REAL = ROOT / "tests" / "fixtures" / "real_cp2k_2026.2"


@pytest.fixture(scope="module")
def rks(tmp_path_factory):
    folder = tmp_path_factory.mktemp("rks101")
    folder, *_ = generate_demo(folder, mesh=(3, 3, 3))
    return folder, load_system(folder / "cohp_mo-1_0.mokp", folder)


def test_fortran_bare_exponent(rks, tmp_path):
    """Fortran ES formats drop the 'E' for 3-digit exponents (Fermi-smeared occupations: 2.9-102)."""
    assert _to_floats("1.5 2.91755220991368-102 -3.0+100 4D-2 -1.0E-01 7") == pytest.approx(
        [1.5, 2.91755220991368e-102, -3.0e100, 0.04, -0.1, 7.0])
    lines = (rks[0] / "cohp_mo-1_0.mokp").read_text().split("\n")
    i = next(n for n, ln in enumerate(lines) if ln.startswith("# MO_COEFF_RE")) - 1   # last OCCUPATIONS line
    toks = lines[i].split()
    assert float(toks[-1]) == 0.0
    lines[i] = "    " + "    ".join(toks[:-1] + ["2.91755220991368-102"])
    f = tmp_path / "smeared.mokp"
    f.write_text("\n".join(lines))
    d = read_mokp(f)
    assert d.occ[(0, 0)][-1] == pytest.approx(2.91755220991368e-102)


def test_trapz_compat():
    """np.trapezoid only exists in NumPy >= 2.0; the descriptors must work with NumPy 1.x as well."""
    g = make_grid(-5, 5, 0.01)
    y = broaden(np.array([0.0]), np.array([2.0]), g, 0.3)
    assert trapz(y, g) == pytest.approx(2.0, rel=1e-6)
    assert curve_descriptors(g, y).integral == pytest.approx(2.0, rel=1e-6)


def test_cohp_gauges(rks):
    """E_cov is invariant under H -> H + dS; the raw AO COHP changes by d*COOP; Loewdin ignores the gauge."""
    sysm = rks[1]
    out = {}
    H0 = dict(sysm.H)
    try:
        for shift_ev in (0.0, 7.0):
            sysm.H = {ks: H + (shift_ev / HARTREE_EV) * sysm.S(ks[0]) for ks, H in H0.items()}
            for gauge in ("raw", "ecov", "eref"):
                eng = BondingEngine(sysm, BondingOptions(basis="ao", cohp_gauge=gauge, cohp_shift=1.0 + shift_ev))
                bonds = eng.make_bond(1, 2) + eng.make_bond(2, 3) + eng.make_bond(1, 2, mode="PAIR")
                res = eng.run(bonds, resolve={0})
                out[(shift_ev, gauge)] = [(r.integrals["cohp"]["total"], r.integrals["coop"]["total"]) for r in res]
                assert res[0].resolved_integrals["cohp"][0].sum() == pytest.approx(
                    res[0].integrals["cohp"]["total"], abs=1e-9)
    finally:
        sysm.H = H0
    n = len(out[(0.0, "raw")])
    for i in range(n):
        c0, o0 = out[(0.0, "raw")][i]
        c1, _ = out[(7.0, "raw")][i]
        assert c1 - c0 == pytest.approx(7.0 * o0, abs=1e-8)                                   # raw: gauge dependent
        assert out[(7.0, "ecov")][i][0] == pytest.approx(out[(0.0, "ecov")][i][0], abs=1e-8)  # E_cov: invariant
        assert out[(7.0, "eref")][i][0] == pytest.approx(out[(0.0, "eref")][i][0], abs=1e-8)  # follows the reference
        assert out[(0.0, "eref")][i][0] == pytest.approx(c0 - 1.0 * o0, abs=1e-8)
    assert BondingEngine(sysm, BondingOptions(basis="lowdin", cohp_gauge="ecov")).gauge == "raw"
    with pytest.raises(ValueError):
        BondingEngine(sysm, BondingOptions(cohp_gauge="nonsense"))


def test_gauge_from_config(tmp_path):
    folder, *_ = generate_demo(tmp_path / "cfg")
    (folder / "cohp.in").write_text("MOKP = cohp_mo-1_0.mokp\nKS_DIR = .\nBOND = 1 2\n")      # ECOV is the default
    sess, reqs, auto = session_from_config(parse_config(folder / "cohp.in"))
    sess.compute(reqs, auto)
    assert sess.engine.gauge == "ecov"
    assert not any(n.startswith("AO 表象的原始 COHP") for n in sess.notes)
    sess.settings.cohp_gauge = "raw"
    sess.compute(reqs, auto)
    assert any(n.startswith("AO 表象的原始 COHP") for n in sess.notes)


def test_fermi_estimate_in_gap(tmp_path):
    """Fermi-Dirac occupations: the estimate must recover mu even when it lies in a gap."""
    folder, *_ = generate_demo(tmp_path / "fd")
    s = BondSession()
    s.load(folder / "cohp_mo-1_0.mokp", folder)
    sysm = s.system
    vbm, cbm, _ = sysm.band_edges()
    mu, kt = vbm + 0.2 * (cbm - vbm), 0.05
    for ks in sysm.d.occ:
        sysm.d.occ[ks] = sysm.g / (1.0 + np.exp((sysm.eps_ev(*ks) - mu) / kt))
    assert sysm.fermi_estimate() == pytest.approx(mu, abs=1e-6)


def test_time_reversal_unfolding(rks):
    """A mesh reduced by time reversal only (CP2K default) is unfolded and reproduces the full mesh."""
    import dataclasses
    from bonddos.kspace import KSystem
    d = rks[1].d
    keep, wts = [], []
    for i, k in enumerate(d.kpoints):
        s = [k + d.kpoints[j] for j in keep]
        if any(np.allclose(x - np.round(x), 0, atol=1e-8) for x in s):
            continue                                   # the −k partner is already kept
        self_inverse = np.allclose(2 * k - np.round(2 * k), 0, atol=1e-8)
        keep.append(i)
        wts.append(d.weights[i] * (1 if self_inverse else 2))
    assert len(keep) == 14
    red = dataclasses.replace(
        d, nkp=len(keep), kpoints=d.kpoints[keep], weights=np.array(wts),
        eig={(n, s): d.eig[(i, s)] for n, i in enumerate(keep) for s in range(d.nspins)},
        occ={(n, s): d.occ[(i, s)] for n, i in enumerate(keep) for s in range(d.nspins)},
        coeff={(n, s): d.coeff[(i, s)] for n, i in enumerate(keep) for s in range(d.nspins)},
        overlap={n: d.overlap[i] for n, i in enumerate(keep)}, warnings=[])
    out = []
    for data in (d, red):
        sysm = KSystem(data)                           # H rebuilt from S C ε C† S when NMO = NAO
        assert sysm.mesh.full and sysm.nk == 27
        eng = BondingEngine(sysm, BondingOptions())
        res = eng.run(eng.make_bond(1, 2, all_images=True) + eng.make_bond(2, 3))
        out.append([r.integrals[q]["total"] for r in res for q in eng.quantities])
    assert any("时间反演" in w for w in sysm.warnings)
    assert out[1] == pytest.approx(out[0], abs=1e-9)


def test_default_gauge_is_ecov():
    from bonddos.session import SessionSettings
    assert SessionSettings().cohp_gauge == "ecov" and BondingOptions().cohp_gauge == "ecov"


def test_cp2k_template():
    """The printed CP2K fragment: required keywords present, and no inline comments (CP2K would take a
    '! comment' after FILENAME as part of the file name)."""
    import re
    from bonddos.cp2k_template import cp2k_template, PRINT_BLOCK
    inc = (ROOT / "examples" / "cp2k_input" / "bonddos_print.inc").read_text(encoding="utf-8")
    assert inc.replace("\r\n", "\n") == PRINT_BLOCK
    t = cp2k_template(mesh=(4, 4, 4), uks=True, smear=True)
    for key in ("&KPOINTS", "GAMMA_CENTERED T", "ADDED_MOS", "&MO_KP ON", "AO_EXPORT_TYPE OVERLAP_MATRIX",
                "&KS_CSR_WRITE ON", "UKS T", "&SMEAR ON"):
        assert key in t
    assert not re.search(r"^\s*[^!\s][^!\n]*!", t, re.M)
    g = cp2k_template(mesh=(1, 1, 1), dos=False)
    assert "SCHEME GAMMA" in g and "GAMMA_CENTERED" not in g and "&DOS" not in g
    assert "GAMMA_CENTERED" not in cp2k_template(mesh=(3, 3, 1))


def test_csr_upper_triangle_autodetect(rks, tmp_path):
    from bonddos.mokp import read_csr, find_ks_sets
    sysm = rks[1]
    f = next(iter(next(iter(find_ks_sets(rks[0]).values())).values()))
    full = read_csr(f, sysm.nao)
    rows = [ln.split() for ln in f.read_text().splitlines() if ln.strip()]
    up = tmp_path / "upper.csr"
    up.write_text("\n".join(" ".join(r) for r in rows if int(r[0]) <= int(r[1])) + "\n")
    assert np.allclose(read_csr(up, sysm.nao), full, atol=1e-12)          # not declared as UPPER


def test_overlap_rebuilt_or_raw_fallback(rks):
    """AO_EXPORT_TYPE GTO_BASIS (no S(k)): rebuild S from a complete MO set, otherwise fall back to RAW."""
    import dataclasses
    from bonddos.kspace import KSystem
    from bonddos.mokp import find_ks_sets
    d = rks[1].d
    d2 = dataclasses.replace(d, overlap={}, ao_export="GTO_BASIS", warnings=["没有 S(k)"])
    sysm = KSystem(d2, next(iter(find_ks_sets(rks[0]).values())))
    if d.nmo == d.nao:
        assert sysm.has_S and np.allclose(sysm.S(0), d.overlap[0], atol=1e-7)
        assert any("重建" in w for w in sysm.warnings)
    else:
        eng = BondingEngine(sysm, BondingOptions(quantities=("cohp",)))
        assert not sysm.has_S and eng.gauge == "raw" and eng.gauge_fallback


def test_missing_input_hints(tmp_path):
    from bonddos.mokp import diagnose_ks_folder, MOKPError
    (tmp_path / "a").mkdir()
    assert "KS_CSR_WRITE ON" in diagnose_ks_folder(tmp_path / "a")[0]
    (tmp_path / "a" / "x-KS_SPIN_1-1_0.csr").write_text("1 1 0.0\n")
    assert "KPOINTS" in diagnose_ks_folder(tmp_path / "a")[0]
    with pytest.raises(MOKPError, match="KPOINTS"):
        load_system(tmp_path / "a" / "none.mokp", tmp_path / "a")


def test_fixed_moment_common_fermi_zero():
    """CP2K puts the zero of all broadened curves at the higher of the two Fermi energies."""
    from bonddos.dosanalysis import DosProject
    from bonddos.dosfiles import DosData
    e = np.linspace(-1, 1, 5)
    pr = DosProject()
    pr.data = [DosData(path="run/x-ALPHA_k1_curve-1.pdos", kind="Mn", spin="alpha", mode="curve", energy=e,
                       channels={"total": e * 0}, zero="fermi", fermi=12.30, iteration=0),
               DosData(path="run/x-BETA_k1_curve-1.pdos", kind="Mn", spin="beta", mode="curve", energy=e,
                       channels={"total": e * 0}, zero="fermi", fermi=12.86, iteration=0)]
    pr._harmonize_fermi()
    assert [d.fermi for d in pr.data] == [12.86, 12.86] and len(pr.notes) == 1


def test_lobster_layout_read_by_pymatgen(rks, tmp_path):
    """Layout details pymatgen's readers rely on (checked against pymatgen with real exports):
    DOSCAR title token 17 = E_Fermi, GROSSPOP has exactly three header lines."""
    s = BondSession()
    s.load(rks[0] / "cohp_mo-1_0.mokp", rks[0])
    s.compute([BondRequest.parse("1 2")])
    s.export(tmp_path, "T")
    dos = (tmp_path / "DOSCAR.lobster").read_text().split("\n")
    assert float(dos[4].split()[17]) == 0.0
    assert int(dos[5].split()[2]) == len(s.grid())
    gp = (tmp_path / "GROSSPOP.lobster").read_text().split("\n")
    first = gp[3].split()
    assert first[0] == "1" and len(first) == 5 and all(ln.startswith("#") for ln in gp[:3])


def test_cli_pdos_cohp_plots_headless(tmp_path, monkeypatch):
    """Everything the GUI does is reachable from the command line (Linux servers without a display)."""
    from bonddos import cli
    folder, *_ = generate_demo(tmp_path / "cli")
    out = tmp_path / "o"
    rc = cli.main(["pdos", str(folder / "cohp_mo-1_0.mokp"), "--resolution", "l", "--out", str(out / "p.csv"),
                   "--stats", str(out / "s.csv"), "--charges", str(out / "q.csv"), "--plot", str(out / "p.png")])
    assert rc == 0 and all((out / n).stat().st_size > 0 for n in ("p.csv", "s.csv", "q.csv", "p.png"))
    (folder / "c.in").write_text("MOKP = cohp_mo-1_0.mokp\nKS_DIR = .\nBOND = 1 2\nOUTPUT_PREFIX = r/X\n")
    assert cli.main(["cohp", str(folder / "c.in"), "--plot"]) == 0
    assert (folder / "r" / "X_COHP.png").stat().st_size > 0 and (folder / "r" / "ICOHPLIST.lobster").exists()
    monkeypatch.setitem(sys.modules, "PySide6", None)          # no GUI toolkit installed
    monkeypatch.delitem(sys.modules, "bonddos.gui.app", raising=False)
    assert cli.main([]) == 0                                   # prints the help instead of crashing


def test_ks_set_selected_by_eigen_equation(tmp_path):
    """Several KS file sets in one folder (e.g. one per GEO_OPT step): pick the one matching the .mokp."""
    folder, *_ = generate_demo(tmp_path / "sets")
    for f in list(folder.glob("*KS_SPIN_*-1_0.csr")):
        rows = f.read_text().splitlines()
        first = rows[0].split()
        first[2] = f"{float(first[2]) + 0.5:.14E}"
        (folder / f.name.replace("-1_0.csr", "-1_9.csr")).write_text("\n".join([" ".join(first)] + rows[1:]) + "\n")
    sysm = load_system(folder / "cohp_mo-1_0.mokp", folder)
    assert any("自动选中" in w and "-1_0.csr" in w for w in sysm.warnings)
    assert all(c.ok for c in sysm.validate())


def test_dos_absolute_zero_header_and_scale(tmp_path):
    from bonddos.dosanalysis import DosProject, DosOptions
    from bonddos.dosfiles import read_dos_file
    grid = "\n".join(f"{e:.4f} {np.exp(-(e - 5.2) ** 2):.6f} 0.0" for e in np.arange(0.0, 10.0, 0.05))
    (tmp_path / "x-curve-1.dos").write_text(
        "# DOS at iteration step i = 0\n# E(Fermi) =     0.000000 a.u. =     0.000000 eV\n"
        "# E(HOCO)  =     0.195649 a.u. =     5.323875 eV\n# Energy zero: ABSOLUTE\n"
        "# Gaussian broadening, FWHM = 0.00367493 a.u. = 0.100000 eV\n#  Energy[eV]       Density     Occupation\n" + grid + "\n")
    (tmp_path / "x-k1_curve-1.pdos").write_text(
        "# K-point projected DOS for atomic kind Si at iteration step i = 0\n# E(Fermi) =     0.204841 a.u. =     5.574013 eV\n"
        "# E(HOCO)  =     0.195649 a.u. =     5.323875 eV\n# Energy zero: ABSOLUTE\n"
        "# Gaussian broadening, FWHM = 0.00367493 a.u. = 0.100000 eV\n# Energy[eV]  total  s  p\n"
        + "\n".join(ln + " 0.0" for ln in grid.split("\n")) + "\n")
    assert read_dos_file(tmp_path / "x-curve-1.dos")[0].fermi is None          # CP2K prints 0 here
    pr = DosProject()
    pr.add_files([tmp_path / "x-curve-1.dos", tmp_path / "x-k1_curve-1.pdos"])
    assert all(d.fermi == pytest.approx(5.574013, abs=1e-5) for d in pr.data)  # taken from the sibling file
    s1, _ = pr.build_series(pr.data[:1], DosOptions(channels=("total",)))
    s2, _ = pr.build_series(pr.data[:1], DosOptions(channels=("total",), scale=2.0))
    assert np.allclose(s2[0].y, 2 * s1[0].y) and s1[0].x[0] == pytest.approx(-5.574013, abs=1e-5)


# ------------------------------------------------------------------ real CP2K 2026.2 output
# tests/fixtures/real_cp2k_2026.2/<case>/ : cohp_mo-1_0.mokp + KS csr files written by CP2K 2026.2
# (conda-forge build) and reference.json = PAIR-mode integrals of the equivalent Γ-point supercell
# calculation (real-space reference that does not depend on any Fourier/wrapping convention).

def _cases():
    return sorted(p.parent for p in REAL.glob("*/reference.json")) if REAL.is_dir() else []


@pytest.mark.parametrize("case", _cases(), ids=lambda p: p.name)
def test_real_cp2k_against_supercell(case):
    import json
    ref = json.loads((case / "reference.json").read_text())
    for basis in ("ao", "lowdin"):
        s = BondSession()
        s.settings.basis = basis
        s.settings.cohp_gauge = "raw"                  # reference.json holds raw-gauge integrals
        s.load(case / "cohp_mo-1_0.mokp", case)
        assert all(c.ok for c in s.system.checks)
        assert s.gauge["available"] and s.gauge["consistent"]
        assert s.eref == pytest.approx(ref["eref_eV"], abs=2e-4)
        res = s.compute([BondRequest(b["a"], b["b"], tuple(b["R"])) for b in ref["bonds"]])
        for b, r in zip(ref["bonds"], res):
            assert r.status == b["status"]
            assert r.bond.length == pytest.approx(b["length"], abs=1e-5)
            for q, vals in b[basis].items():
                for key, v in vals.items():
                    assert r.integrals[q][key] == pytest.approx(v, abs=ref["tol"]), (basis, b, q, key)
        if basis == "ao":
            pm, _ = s.compute_populations()
            assert pm.charges == pytest.approx(ref["mulliken_charges_cp2k"], abs=2e-5)
            if "mulliken_moments_cp2k" in ref:
                assert pm.moments == pytest.approx(ref["mulliken_moments_cp2k"], abs=2e-5)
