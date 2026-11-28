"""Self tests: python -m pytest tests -q   (run from the project folder)."""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

from bonddos.demo import generate_demo, Model, mp_mesh, solve, write_mokp, write_csr  # noqa: E402
from bonddos.kspace import load_system, infer_mesh  # noqa: E402
from bonddos.bonding import BondingEngine, BondingOptions, gauge_check  # noqa: E402
from bonddos.session import BondSession, BondRequest, parse_config, session_from_config  # noqa: E402
from bonddos.population import populations  # noqa: E402
from bonddos.dosfiles import read_dos_file  # noqa: E402
from bonddos.dosanalysis import DosProject, DosOptions  # noqa: E402
from bonddos.spectra import broaden, make_grid, curve_descriptors, trapz  # noqa: E402
from supercell_ref import reference_integrals  # noqa: E402


@pytest.fixture(scope="module")
def rks(tmp_path_factory):
    folder = tmp_path_factory.mktemp("rks")
    folder, model, kpts, ef = generate_demo(folder, mesh=(3, 3, 3))
    return folder, model, load_system(folder / "cohp_mo-1_0.mokp", folder)


def test_checks_pass(rks):
    _, _, sysm = rks
    for c in sysm.validate():
        assert c.ok, c


@pytest.mark.parametrize("basis", ["ao", "lowdin"])
def test_against_supercell(rks, basis):
    _, model, sysm = rks
    eng = BondingEngine(sysm, BondingOptions(basis=basis, cohp_gauge="raw"))
    bonds = []
    for a, b in [(1, 2), (1, 3), (2, 3), (1, 1)]:
        bonds += eng.make_bond(a, b, all_images=True)
    for r in eng.run(bonds):
        ref = reference_integrals(model, (3, 3, 3), r.bond.a, r.bond.b, r.bond.r_wrap)
        if basis == "ao":
            assert r.integrals["cohp"]["total"] == pytest.approx(ref["cohp"], rel=1e-6, abs=1e-7)
            assert r.integrals["coop"]["total"] == pytest.approx(ref["coop"], rel=1e-6, abs=1e-8)
            assert r.integrals["cobi"]["total"] == pytest.approx(ref["cobi_m"], rel=1e-6, abs=1e-8)
        else:
            assert r.integrals["cohp"]["total"] == pytest.approx(ref["cohp_l"], rel=1e-6, abs=1e-7)
            assert r.integrals["cobi"]["total"] == pytest.approx(ref["cobi_l"], rel=1e-6, abs=1e-8)


def test_orbital_resolution_sums(rks):
    _, _, sysm = rks
    for res in ("l", "lm", "shell", "ao"):
        eng = BondingEngine(sysm, BondingOptions(basis="ao", resolution=res))
        bonds = eng.make_bond(1, 2)
        r = eng.run(bonds, resolve={0})[0]
        for q in ("cohp", "coop", "cobi"):
            assert r.resolved_integrals[q][0].sum() == pytest.approx(r.integrals[q]["total"], abs=1e-9)


def test_pair_equals_sum_of_images(rks):
    _, _, sysm = rks
    eng = BondingEngine(sysm, BondingOptions(basis="ao", quantities=("cohp", "coop")))
    st = sysm.structure
    imgs = [eng.make_bond(1, 2, R)[0] for R, d in st.images(1, 2, 12.0)]
    # keep only one representative per alias class of the 3x3x3 mesh
    seen, uniq = set(), []
    for b in imgs:
        key = tuple(x % 3 for x in b.r_wrap)
        if key not in seen:
            seen.add(key)
            uniq.append(b)
    pair = eng.make_bond(1, 2, mode="PAIR")
    res = eng.run(uniq + pair)
    tot = sum(r.integrals["cohp"]["total"] for r in res[:-1])
    assert tot == pytest.approx(res[-1].integrals["cohp"]["total"], abs=1e-8)


def test_gauge_invariance(tmp_path):
    """Raw coordinates shifted by lattice vectors must give identical bonds (with relabelled R)."""
    out = []
    for shift in (False, True):
        m = Model(raw_shift=shift)
        k = mp_mesh((3, 3, 3))
        eig, C, occ, _ = solve(m, k)
        f = tmp_path / f"s{int(shift)}"
        f.mkdir()
        write_mokp(f / "a.mokp", m, k, eig, C, occ)
        write_csr(f, m, k)
        sysm = load_system(f / "a.mokp", f)
        eng = BondingEngine(sysm, BondingOptions())
        bonds = eng.make_bond(2, 3) + eng.make_bond(1, 3)
        out.append([(round(b.length, 6), r.integrals["cohp"]["total"]) for b, r in zip(bonds, eng.run(bonds))])
    for (d0, v0), (d1, v1) in zip(*out):
        assert d0 == d1
        assert v0 == pytest.approx(v1, abs=1e-9)


def test_gauge_check(rks):
    g = gauge_check(rks[2])
    assert g["available"] and g["consistent"]


def test_populations(rks):
    _, model, sysm = rks
    for method in ("mulliken", "lowdin"):
        p = populations(sysm, method)
        assert p.atom_pop.sum() == pytest.approx(model.nelectrons(), abs=1e-6)


def test_uks_and_shifted_mesh(tmp_path):
    folder, model, kpts, _ = generate_demo(tmp_path / "uks", uks=True, mesh=(4, 4, 4), gamma_centered=False)
    s = BondSession()
    s.load(folder / "cohp_mo-1_0.mokp", folder)
    assert s.system.mesh.mesh == (4, 4, 4) and s.system.mesh.full
    assert all(c.ok for c in s.system.checks)
    res = s.compute([BondRequest.parse("1 2 ORB")], auto={"dmin": 1.0, "dmax": 2.2})
    r = res[0]
    assert set(r.integrals["cohp"]) == {"up", "down", "total"}
    assert r.integrals["cohp"]["total"] == pytest.approx(r.integrals["cohp"]["up"] + r.integrals["cohp"]["down"])
    files = s.export(tmp_path / "out", "X")
    assert any(p.name == "ICOHPLIST.lobster" for p in files)


def test_spectrum_integral_consistency(rks):
    """For an insulator the cumulative (capacity) curve at E_F equals the occupied integral."""
    folder, _, _ = rks
    s = BondSession()
    s.load(folder / "cohp_mo-1_0.mokp", folder)
    gap = s.summary()["gap"]
    assert gap > 0.05
    s.settings.emin, s.settings.emax, s.settings.de = -40, 0.5 * gap, 0.002
    s.settings.fwhm = 0.02
    r = s.compute([BondRequest.parse("1 2")])[0]
    g, ci = s.curves(r, "cohp", cumulative=True)
    assert ci[0][-1] == pytest.approx(r.integrals["cohp"]["total"], abs=2e-3)


def test_config_roundtrip(tmp_path):
    folder, *_ = generate_demo(tmp_path / "cfg")
    (folder / "cohp.in").write_text("MOKP = cohp_mo-1_0.mokp\nKS_DIR = .\nSIGMA_EV = 0.1\nBOND = 1 2\n"
                                    "BOND = 2 3 ALL\nBOND = 1 2 PAIR\nAUTO_BONDS = 1.0 2.0\nEF_EV = auto\n")
    sess, reqs, auto = session_from_config(parse_config(folder / "cohp.in"))
    assert len(reqs) == 3 and auto["dmax"] == 2.0
    assert sess.settings.fwhm == pytest.approx(0.235482, rel=1e-5)


def test_lobster_files_parse(tmp_path):
    folder, *_ = generate_demo(tmp_path / "lob")
    s = BondSession()
    s.load(folder / "cohp_mo-1_0.mokp", folder)
    s.compute([BondRequest.parse("1 2"), BondRequest.parse("2 3")])
    s.export(tmp_path / "o", "T")
    lines = (tmp_path / "o" / "COHPCAR.lobster").read_text().split("\n")
    par = lines[1].split()
    nb = int(par[0]) - 1
    assert nb == 2 and int(par[1]) == 1
    data = np.array([row.split() for row in lines[nb + 3:] if row.strip()], float)
    assert data.shape[1] == 1 + 2 * (nb + 1)
    ic = [ln.split() for ln in (tmp_path / "o" / "ICOHPLIST.lobster").read_text().split("\n")[1:-1]]
    assert all(len(t) == 8 for t in ic)
    ch = (tmp_path / "o" / "CHARGE.lobster").read_text().split("\n")[3:-3]
    assert len(ch) == 3 and all(len(t.split()) == 4 for t in ch)


def test_mesh_inference():
    k = mp_mesh((2, 3, 1), gamma_centered=False)
    info = infer_mesh(k, np.full(len(k), 1 / len(k)))
    assert info.mesh == (2, 3, 1) and info.full


def test_dos_formats():
    fx = ROOT / "tests" / "fixtures"
    a = read_dos_file(fx / "Fe-ALPHA_k1_curve-1.pdos")[0]
    assert a.kpoint and a.spin == "alpha" and a.zero == "fermi" and a.fwhm == pytest.approx(0.1)
    assert trapz(a.channels["d"], a.energy) == pytest.approx(5.0, abs=1e-3)
    t = read_dos_file(fx / "FeO-curve-1.dos")
    assert [x.spin for x in t] == ["alpha", "beta"]
    r = read_dos_file(fx / "TEST-k2-1.pdos")[0]
    assert r.mode == "sticks" and "p" in r.channels
    assert np.allclose(r.channels["p"], r.channels["py"] + r.channels["pz"] + r.channels["px"])


def test_wo3_band_center_and_merge():
    pr = DosProject()
    pr.add_files(sorted((ROOT / "examples" / "WO3_gamma_pdos").glob("*.pdos")))
    assert not pr.errors and len(pr.data) == 2
    series, _ = pr.build_series(pr.data, DosOptions(channels=("total", "d"), fwhm=0.1))
    rows = pr.statistics(series)
    w_total = [r for r in rows if r["name"] == "W" and r["channel"] == "total"][0]
    # discrete statistics use raw MO weights (independent of broadening)
    e, wt = series[0].sticks
    assert w_total["integral"] == pytest.approx(wt.sum())
    merged, _ = pr.build_series(pr.data, DosOptions(channels=("total",), merge=True))
    assert len(merged) == 1
    assert merged[0].sticks[1].sum() == pytest.approx(2500.0, rel=1e-6)  # one state per MO


def test_broaden_area():
    g = make_grid(-10, 10, 0.005)
    y = broaden(np.array([0.0, 1.0]), np.array([1.0, 2.0]), g, 0.2)
    assert trapz(y, g) == pytest.approx(3.0, rel=1e-6)
    d = curve_descriptors(g, y)
    assert d.center == pytest.approx(2.0 / 3.0, abs=1e-4)
