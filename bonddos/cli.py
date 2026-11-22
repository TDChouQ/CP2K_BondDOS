"""Command line interface.

  CP2K_BondDOS                         -> start the GUI
  CP2K_BondDOS cohp cohp.in [--check]  -> k-point COHP/COOP/COBI, charges, PDOS, LOBSTER files
  CP2K_BondDOS dos *.pdos [options]    -> CP2K DOS/PDOS curves and band statistics
  CP2K_BondDOS merge_pdos *.pdos       -> merged_pdos.dat
  CP2K_BondDOS band_center *.pdos      -> band_center.dat
  CP2K_BondDOS demo <folder>           -> synthetic MO_KP/KS test data
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__


def _progress(f, text):
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
    for w in info["warnings"]:
        print("  警告:", w)
    if args.check:
        return 0 if all(c.ok for c in sess.system.checks) else 2
    res = sess.compute(reqs, auto, progress=_progress)
    sys.stderr.write("\n")
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
    opt = DosOptions(reference=args.ref, fwhm=args.fwhm, channels=chans, merge=args.merge, occupied=args.occupied)
    series, warns = pr.build_series(pr.data, opt)
    for w in warns:
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
    sub = p.add_subparsers(dest="cmd")
    c = sub.add_parser("cohp", help="k 点 COHP/COOP/COBI 分析")
    c.add_argument("config")
    c.add_argument("--check", action="store_true", help="只检查数据")
    c.set_defaults(func=cmd_cohp)
    d = sub.add_parser("dos", help="CP2K DOS/PDOS 曲线与能带中心")
    d.add_argument("files", nargs="+")
    d.add_argument("--ref", default="fermi", choices=["file", "fermi", "hoco", "absolute"])
    d.add_argument("--fwhm", type=float, default=0.1)
    d.add_argument("--channels", default="total,s,p,d,f")
    d.add_argument("--merge", action="store_true")
    d.add_argument("--occupied", action="store_true")
    d.add_argument("--emin", type=float)
    d.add_argument("--emax", type=float)
    d.add_argument("--out")
    d.add_argument("--stats")
    d.add_argument("--plot")
    d.set_defaults(func=cmd_dos)
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
    e = sub.add_parser("demo", help="生成合成演示数据")
    e.add_argument("folder")
    e.add_argument("--uks", action="store_true")
    e.add_argument("--metal", action="store_true")
    e.set_defaults(func=cmd_demo)
    return p


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        from .gui.app import run_gui
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
