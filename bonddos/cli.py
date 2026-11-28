"""Command line interface.

  CP2K_BondDOS                         -> start the GUI
  CP2K_BondDOS cohp cohp.in [--check]  -> k-point COHP/COOP/COBI, charges, PDOS, LOBSTER files
  CP2K_BondDOS pdos x.mokp --group "Mn: Mn" --group "O: 2,5-7"  -> atom/orbital PDOS, band centres, charges
  CP2K_BondDOS dos *.pdos [options]    -> CP2K DOS/PDOS curves and band statistics
  CP2K_BondDOS merge_pdos *.pdos       -> merged_pdos.dat
  CP2K_BondDOS band_center *.pdos      -> band_center.dat
  CP2K_BondDOS demo <folder>           -> synthetic MO_KP/KS test data
  CP2K_BondDOS template [--uks --smear --mesh 3 3 3]  -> what to add to the CP2K input
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__


def _progress(f, text):
    if not sys.stderr.isatty():      # redirected to a log file / job script: no carriage-return spam
        return
    sys.stderr.write(f"\r[{f * 100:5.1f}%] {text[:60]:60s}")
    sys.stderr.flush()


def cmd_cohp(args):
    from .session import parse_config, session_from_config
    cfg = parse_config(args.config)
    sess, reqs, auto = session_from_config(cfg, progress=_progress)
    sys.stderr.write("\n")
    info = sess.summary()
    print(f"体系: {info['formula']}  原子 {info['atoms']}  AO {info['nao']}  MO {info['nmo']}  "
          f"k 点 {info['nkp']} (网格 {info['mesh']})  自旋 {info['nspins']}  H: {info['H']}")
    print(f"能量零点: {sess.eref:.6f} eV ({sess.eref_note})")
    for c in sess.system.checks:
        print(f"  [{'OK' if c.ok else '!!'}] {c.name}: {c.value:.3e} {c.note if not c.ok else ''}")
    if sess.gauge.get("available"):
        print(f"  相位约定检查: score(+)={sess.gauge['score_plus']:.3g} score(−)={sess.gauge['score_minus']:.3g} "
              f"{'一致' if sess.gauge['consistent'] else '不一致!'}")
    for w in info["mesh_messages"] + info["warnings"]:
        print("  警告:", w)
    if args.check:
        return 0 if all(c.ok for c in sess.system.checks) else 2
    res = sess.compute(reqs, auto, progress=_progress)
    sys.stderr.write("\n")
    for w in sess.notes:
        if w not in info["warnings"]:
            print("  警告:", w)
    if sess.engine.gauge != "raw":
        print(f"COHP 规范: {sess.engine.gauge.upper()}")
    for i, r in enumerate(res, 1):
        I = r.integrals
        parts = [f"{k.upper()}={v['total']:.4f}" for k, v in I.items()]
        print(f"{i:4d} {r.bond.label:28s} d={r.bond.length:.4f} Å  {'  '.join(parts)}  {r.status}")
    prefix = cfg.get("OUTPUT_PREFIX", "results/CP2K")
    out = Path(cfg["_dir"]) / prefix
    written = sess.export(out.parent, out.name, lobster_files=cfg.get("LOBSTER_FILES", "T").upper().startswith("T"),
                          progress=_progress)
    sys.stderr.write("\n")
    print(f"已写出 {len(written)} 个文件到 {out.parent}")
    if args.plot:
        for p in _save_cohp_plots(sess, out.parent, out.name):
            print("图片 ->", p)
    return 0


def _legend_outside(fig):
    """Put one de-duplicated legend to the right of the first panel (keeps the curves unobstructed)."""
    ax = fig.axes[0]
    handles, labels = ax.get_legend_handles_labels()
    for a in fig.axes:
        if a.get_legend() is not None:
            a.get_legend().remove()
    uniq = {}
    for h, lab in zip(handles, labels):
        uniq.setdefault(lab, h)
    if uniq:
        ax.legend(list(uniq.values()), list(uniq.keys()), loc="upper left", bbox_to_anchor=(1.01, 1.0),
                  frameon=False, fontsize=9, borderaxespad=0.0)


def _save_cohp_plots(sess, folder, prefix, max_bonds=8):
    """PNG per quantity (−COHP / COOP / COBI + integrals); symmetry-equivalent bonds are drawn once."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .bonding import QLABEL
    from .plotting import draw_cohp, PlotStyle, apply_rc
    apply_rc()
    st = sess.system.structure
    pick, seen = [], set()
    for r in sess.results:
        if r.status in ("aliased", "nonperiodic"):
            continue
        key = (r.bond.mode, tuple(sorted((st.symbol(r.bond.a), st.symbol(r.bond.b)))), round(r.bond.length, 2))
        if key not in seen:
            seen.add(key)
            pick.append(r)
    spin_name = {0: "up", 1: "down"} if sess.system.ns == 2 else {0: "total"}
    out = []
    for q in sess.engine.quantities:
        curves, cumul = [], []
        for r in pick[:max_bonds]:
            g, c = sess.curves(r, q)
            _, ci = sess.curves(r, q, cumulative=True)
            lab = f"{r.bond.label} ({r.bond.length:.2f} Å)"
            for sp in c:
                curves.append(dict(label=lab, x=g, y=c[sp], spin=spin_name[sp]))
                cumul.append(dict(label=lab, x=g, y=ci[sp], spin=spin_name[sp]))
        if not curves:
            continue
        fig = plt.figure(figsize=(6.4, 5.6), dpi=150)
        draw_cohp(fig, curves, cumul, PlotStyle(xmin=sess.settings.emin, xmax=sess.settings.emax), QLABEL[q], True, True)
        p = Path(folder) / f"{prefix}_{QLABEL[q]}.png"
        _legend_outside(fig)
        fig.savefig(p, dpi=200, bbox_inches="tight")
        plt.close(fig)
        out.append(p)
    return out


def cmd_pdos(args):
    """Atom / orbital resolved k-point PDOS, band descriptors and charges from a .mokp file (no GUI needed)."""
    import csv
    import numpy as np
    from .kspace import parse_atom_selection
    from .population import PDOSGroup, group_ao_channels, stick_weights
    from .session import BondSession
    from .spectra import stick_descriptors
    sess = BondSession()
    s = sess.settings
    s.fwhm, s.emin, s.emax, s.de = args.fwhm, args.emin, args.emax, args.de
    ref = args.ref.lower()
    if ref in ("auto", "hoco", "midgap", "fermi"):
        s.eref_mode = ref
    else:
        s.eref_mode, s.eref_value = "manual", float(ref)
    sess.load(args.mokp, None, periodic=args.periodic, progress=_progress)
    sys.stderr.write("\n")
    sysm, st = sess.system, sess.system.structure
    info = sess.summary()
    print(f"体系: {info['formula']}  AO {info['nao']}  MO {info['nmo']}  k 点 {info['nkp']} (网格 {info['mesh']})  自旋 {info['nspins']}")
    print(f"能量零点: {sess.eref:.6f} eV ({sess.eref_note})")
    for w in info["mesh_messages"] + info["warnings"]:
        if "KS" not in w:            # H(k) is not needed for PDOS / charges
            print("  警告:", w)
    if not sysm.has_S:
        print("错误: 需要 S(k)（MO_KP 的 AO_EXPORT_TYPE OVERLAP_MATRIX）", file=sys.stderr)
        return 1
    groups = []
    specs = args.group or [f"{el}: {el}" for el in dict.fromkeys(a.element for a in sysm.d.atoms)]
    for spec in specs:
        name, _, sel = spec.partition(":")
        sel = sel.strip() or name.strip()
        groups.append(PDOSGroup(name.strip(), parse_atom_selection(sel, len(sysm.d.atoms), st), args.resolution))
    grid, cur = sess.pdos(groups, args.method, True, occupied=args.occupied, progress=_progress)
    sys.stderr.write("\n")
    spin_name = {0: "up", 1: "down"} if sysm.ns == 2 else {0: "total"}
    names = list(cur[0].keys())
    rows = []
    print(f"{'series':24s}{'integral':>10s}{'center':>10s}{'width':>9s}{'occ_center':>12s}{'filling':>9s}   (能窗内，两自旋合计，离散态精确求和)")
    for g in groups:
        for ch, _ in group_ao_channels(sysm, g):
            e, w = stick_weights(sysm, g, args.method, ch)
            d = stick_descriptors(e, w, args.emin, args.emax, zero=sess.eref)
            nm = g.name if ch == "total" else f"{g.name}·{ch}"
            rows.append({"series": nm, "atoms": len(g.atoms), **d.as_dict()})
            print(f"{nm:24s}{d.integral:10.4f}{d.center:10.4f}{d.width:9.4f}{d.occ_center:12.4f}{d.filling:9.4f}")
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8-sig") as fh:
        wr = csv.writer(fh)
        wr.writerow(["E-Eref_eV"] + [f"{n}[{spin_name[sp]}]" for sp in cur for n in names])
        for row in np.column_stack([grid] + [cur[sp][n] for sp in cur for n in names]):
            wr.writerow([f"{x:.6g}" for x in row])
    print("PDOS 曲线 ->", out)
    if args.stats:
        keys = ["series", "atoms", "integral", "center", "width", "skewness", "kurtosis", "filling", "occ_integral", "occ_center"]
        with open(args.stats, "w", newline="", encoding="utf-8-sig") as fh:
            wr = csv.writer(fh)
            wr.writerow(keys)
            for r in rows:
                wr.writerow([r[k] if isinstance(r[k], (str, int)) else f"{r[k]:.6f}" for k in keys])
        print("能带统计 ->", args.stats)
    if args.charges:
        pm, pl = sess.compute_populations()
        with open(args.charges, "w", newline="", encoding="utf-8-sig") as fh:
            wr = csv.writer(fh)
            wr.writerow(["atom", "element", "Z_eff", "Mulliken_charge", "Loewdin_charge"] +
                        (["Mulliken_moment", "Loewdin_moment"] if sysm.ns == 2 else []))
            for i, a in enumerate(sysm.d.atoms):
                wr.writerow([a.index, a.element, a.z_eff, f"{pm.charges[i]:.5f}", f"{pl.charges[i]:.5f}"] +
                            ([f"{pm.moments[i]:.5f}", f"{pl.moments[i]:.5f}"] if sysm.ns == 2 else []))
        print("电荷 ->", args.charges)
    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from .plotting import draw_lines, PlotStyle, apply_rc
        apply_rc()
        spin_plot = {0: "alpha", 1: "beta"} if sysm.ns == 2 else {0: "none"}
        series = [dict(label=n, x=grid, y=cur[sp][n], spin=spin_plot[sp]) for sp in cur for n in names
                  if n != "TDOS" or args.tdos]
        fig = plt.figure(figsize=(6.4, 4.2), dpi=150)
        draw_lines(fig, series, PlotStyle(xmin=args.emin, xmax=args.emax))
        _legend_outside(fig)
        fig.savefig(args.plot, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print("图片 ->", args.plot)
    return 0


def cmd_dos(args):
    import numpy as np
    from .dosanalysis import DosProject, DosOptions, write_series_csv, write_stats_csv
    pr = DosProject()
    pr.add_files(args.files)
    for e in pr.errors:
        print("错误:", e)
    if not pr.data:
        return 1
    chans = tuple(x.strip() for x in args.channels.split(","))
    opt = DosOptions(reference=args.ref, fwhm=args.fwhm, channels=chans, merge=args.merge, occupied=args.occupied,
                     scale=args.scale)
    series, warns = pr.build_series(pr.data, opt)
    for w in pr.notes + warns:
        print("警告:", w)
    rows = pr.statistics(series, args.emin, args.emax)
    print(f"{'name':20s}{'ch':>8s}{'spin':>7s}{'integral':>12s}{'center':>11s}{'width':>9s}{'occ_center':>12s}{'filling':>9s}")
    for r in rows:
        print(f"{r['name']:20s}{r['channel']:>8s}{r['spin']:>7s}{r['integral']:12.4f}{r['center']:11.4f}"
              f"{r['width']:9.4f}{r['occ_center']:12.4f}{r['filling']:9.4f}")
    if args.out:
        write_series_csv(args.out, series)
        print("曲线 ->", args.out)
    if args.stats:
        write_stats_csv(args.stats, rows)
        print("统计 ->", args.stats)
    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from .plotting import draw_lines, PlotStyle, apply_rc
        apply_rc()
        fig = plt.figure(figsize=(6.4, 4.2), dpi=150)
        draw_lines(fig, [dict(label=f"{s.name}·{s.channel}", x=s.x, y=s.y, spin=s.spin) for s in series],
                   PlotStyle(xmin=args.emin, xmax=args.emax))
        fig.savefig(args.plot, dpi=300)
        print("图片 ->", args.plot)
    return 0


def cmd_merge(args):
    from .dosanalysis import DosProject, DosOptions, write_series_csv
    pr = DosProject()
    pr.add_files(args.files)
    for e in pr.errors:
        print("错误:", e)
    chans = tuple(c for c in pr.available_channels() if c in ("total", "s", "p", "d", "f"))
    series, warns = pr.build_series(pr.data, DosOptions(reference=args.ref, fwhm=args.fwhm, channels=chans, merge=True))
    write_series_csv(args.out, series)
    print(f"合并 {len(pr.data)} 组数据 -> {args.out}")
    return 0


def cmd_band_center(args):
    from .dosanalysis import DosProject, DosOptions
    pr = DosProject()
    pr.add_files(args.files)
    chans = tuple(c for c in pr.available_channels() if c in ("total", "s", "p", "d", "f"))
    series, _ = pr.build_series(pr.data, DosOptions(reference=args.ref, fwhm=args.fwhm, channels=chans,
                                                    merge=len(pr.data) > 1))
    rows = pr.statistics(series, args.emin, args.emax)
    lines = ["=" * 62, f"{'Band Center (eV)':^62s}", "=" * 62,
             f"{'Orbital':10s}{'Spin':8s}{'All':>14s}{'Occupied(E<=0)':>18s}{'Width':>10s}"]
    for r in rows:
        lines.append(f"{r['channel']:10s}{r['spin']:8s}{r['center']:14.6f}{r['occ_center']:18.6f}{r['width']:10.4f}")
    text = "\n".join(lines)
    print(text)
    Path(args.out).write_text(text + "\n", encoding="utf-8")
    return 0


def cmd_template(args):
    from .cp2k_template import cp2k_template
    text = cp2k_template(mesh=tuple(args.mesh), added_mos=args.added_mos, uks=args.uks, smear=args.smear,
                         dos=not args.no_dos)
    if args.out:
        Path(args.out).write_text(text, encoding="utf-8")
        print("CP2K 输入片段 ->", args.out)
    else:
        sys.stdout.write(text)
    return 0


def cmd_demo(args):
    from .demo import generate_demo
    folder, *_ = generate_demo(args.folder, uks=args.uks, metal=args.metal)
    cfg = Path(folder) / "cohp.in"
    cfg.write_text("MOKP = cohp_mo-1_0.mokp\nKS_DIR = .\nPERIODIC = XYZ\nCOHP_BASIS = AO\nORBITAL_RESOLUTION = L\n"
                   "FWHM_EV = 0.2\nEMIN_EV = -20\nEMAX_EV = 10\nEF_EV = auto\nBOND = 1 2 ORB\nBOND = 2 3 ALL\n"
                   "AUTO_BONDS = 1.0 2.6\nOUTPUT_PREFIX = results/DEMO\n", encoding="utf-8")
    print("演示数据 ->", folder)
    return 0


def build_parser():
    p = argparse.ArgumentParser(prog="CP2K_BondDOS", description=f"CP2K BondDOS {__version__}")
    p.add_argument("--version", action="version", version=f"CP2K BondDOS {__version__}")
    sub = p.add_subparsers(dest="cmd")
    c = sub.add_parser("cohp", help="k 点 COHP/COOP/COBI 分析")
    c.add_argument("config")
    c.add_argument("--check", action="store_true", help="只检查数据")
    c.add_argument("--plot", action="store_true", help="同时保存 −COHP/COOP/COBI 曲线图（PNG，对称等价的键只画一根）")
    c.set_defaults(func=cmd_cohp)
    d = sub.add_parser("dos", help="CP2K DOS/PDOS 曲线与能带中心")
    d.add_argument("files", nargs="+")
    d.add_argument("--ref", default="fermi", choices=["file", "fermi", "hoco", "absolute"])
    d.add_argument("--fwhm", type=float, default=0.1)
    d.add_argument("--channels", default="total,s,p,d,f")
    d.add_argument("--merge", action="store_true")
    d.add_argument("--occupied", action="store_true")
    d.add_argument("--scale", type=float, default=1.0,
                   help="强度乘以该系数；CP2K 闭壳层(RKS)曲线每个轨道计 1 个态，要含两个自旋时用 2")
    d.add_argument("--emin", type=float)
    d.add_argument("--emax", type=float)
    d.add_argument("--out")
    d.add_argument("--stats")
    d.add_argument("--plot")
    d.set_defaults(func=cmd_dos)
    a = sub.add_parser("pdos", help="由 .mokp 计算任意原子/原子组的 PDOS、能带中心和电荷（不需拆分 KIND）")
    a.add_argument("mokp")
    a.add_argument("--group", action="append", metavar="'名称: 选择'",
                   help="原子组，可重复；选择可写元素符号或编号区间，如 'Mn_surf: 1-4,9'。不给时每种元素一组")
    a.add_argument("--resolution", default="l", choices=["total", "l", "lm", "shell"])
    a.add_argument("--method", default="mulliken", choices=["mulliken", "lowdin"])
    a.add_argument("--fwhm", type=float, default=0.1)
    a.add_argument("--emin", type=float, default=-15.0)
    a.add_argument("--emax", type=float, default=10.0)
    a.add_argument("--de", type=float, default=0.01)
    a.add_argument("--ref", default="auto", help="auto | hoco | midgap | fermi | 绝对能量(eV)")
    a.add_argument("--periodic", default="XYZ")
    a.add_argument("--occupied", action="store_true", help="按占据数加权（电子/eV）")
    a.add_argument("--tdos", action="store_true", help="图中同时画总 DOS")
    a.add_argument("--out", default="atom_pdos.csv")
    a.add_argument("--stats")
    a.add_argument("--charges")
    a.add_argument("--plot")
    a.set_defaults(func=cmd_pdos)
    m = sub.add_parser("merge_pdos", help="合并 PDOS -> merged_pdos.dat")
    m.add_argument("files", nargs="+")
    m.add_argument("--ref", default="fermi")
    m.add_argument("--fwhm", type=float, default=0.1)
    m.add_argument("--out", default="merged_pdos.dat")
    m.set_defaults(func=cmd_merge)
    b = sub.add_parser("band_center", help="能带中心 -> band_center.dat")
    b.add_argument("files", nargs="+")
    b.add_argument("--ref", default="fermi")
    b.add_argument("--fwhm", type=float, default=0.1)
    b.add_argument("--emin", type=float)
    b.add_argument("--emax", type=float)
    b.add_argument("--out", default="band_center.dat")
    b.set_defaults(func=cmd_band_center)
    t = sub.add_parser("template", help="打印 CP2K 输入中需要加入的内容（&KPOINTS / ADDED_MOS / &MO_KP / &KS_CSR_WRITE）")
    t.add_argument("--mesh", type=int, nargs=3, default=[3, 3, 3], metavar=("N1", "N2", "N3"))
    t.add_argument("--added-mos", type=int, default=50)
    t.add_argument("--uks", action="store_true")
    t.add_argument("--smear", action="store_true")
    t.add_argument("--no-dos", action="store_true", help="不包含可选的 &DOS 段")
    t.add_argument("--out")
    t.set_defaults(func=cmd_template)
    e = sub.add_parser("demo", help="生成合成演示数据")
    e.add_argument("folder")
    e.add_argument("--uks", action="store_true")
    e.add_argument("--metal", action="store_true")
    e.set_defaults(func=cmd_demo)
    return p


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        try:
            from .gui.app import run_gui
        except ImportError:        # headless install (Linux server): no PySide6
            build_parser().print_help()
            print("\n没有安装图形界面组件（PySide6）。无界面环境请直接使用上面的子命令；"
                  "需要界面时 pip install PySide6-Essentials")
            return 0
        return run_gui()
    args = build_parser().parse_args(argv)
    if not getattr(args, "func", None):
        build_parser().print_help()
        return 1
    try:
        return args.func(args)
    except Exception as e:  # pragma: no cover
        print(f"\n错误: {e}", file=sys.stderr)
        return 1
