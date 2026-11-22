"""Turn CP2K .dos/.pdos data into plot series and band statistics."""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .dosfiles import DosData, read_dos_file
from .spectra import broaden, make_grid, curve_descriptors, stick_descriptors, cumulative_trapz


@dataclass(eq=False)
class Series:
    name: str
    spin: str
    x: np.ndarray
    y: np.ndarray
    source: str = ""
    channel: str = "total"
    kind: str = ""
    sticks: tuple | None = None      # (energies, weights) when derived from raw per-MO data
    note: str = ""


@dataclass
class DosOptions:
    reference: str = "fermi"      # file|fermi|hoco|absolute
    shift: float = 0.0            # extra shift added to energies (eV)
    fwhm: float = 0.10            # re-broadening of raw (stick) data, eV
    shape: str = "gaussian"
    mixing: float = 0.5
    de: float = 0.01
    channels: tuple = ("total",)
    merge: bool = False           # sum selected datasets per spin & channel
    occupied: bool = False        # occupation-weighted (raw data only)
    rebroaden_curves: bool = False
    integrate: bool = False       # IDOS instead of DOS


class DosProject:
    def __init__(self):
        self.data: list[DosData] = []
        self.errors: list[str] = []

    def add_files(self, paths):
        new = []
        for p in paths:
            try:
                ds = read_dos_file(p)
            except Exception as e:  # report and continue
                self.errors.append(f"{Path(p).name}: {e}")
                continue
            for d in ds:
                if any(Path(x.path) == Path(d.path) and x.spin == d.spin and x.kind == d.kind for x in self.data):
                    continue
                self.data.append(d)
                new.append(d)
        return new

    def clear(self):
        self.data.clear()
        self.errors.clear()

    def available_channels(self, items=None):
        items = self.data if items is None else items
        seen = []
        for d in items:
            for k in d.channels:
                if k not in seen and k != "occupied":
                    seen.append(k)
        order = ["total", "s", "p", "d", "f", "g"]
        return [c for c in order if c in seen] + [c for c in seen if c not in order]

    # ------------------------------------------------------------------ series
    def build_series(self, items: list[DosData], opt: DosOptions):
        out = []
        warn = []
        for d in items:
            e = d.energies_relative_to(opt.reference)
            if e is None:
                warn.append(f"{d.name}: 文件中没有 {opt.reference.upper()} 参考能级，按文件原能量零点显示")
                e = d.energy.copy()
            e = e + opt.shift
            for ch in opt.channels:
                if ch not in d.channels:
                    continue
                y = d.channels[ch]
                if d.mode == "sticks":
                    # one state per MO; occupation-weighted spectrum counts electrons (f <= 2 or 1)
                    w = y * d.occupation if opt.occupied else y
                    fw = max(opt.fwhm, 0.0)
                    lo = e.min() - 6 * max(fw, opt.de)
                    hi = e.max() + 6 * max(fw, opt.de)
                    grid = make_grid(lo, hi, opt.de)
                    curve = broaden(e, w, grid, fw if fw > 0 else 0.0, opt.shape, opt.mixing)
                    out.append(Series(d.label, d.spin, grid, curve, d.name, ch, d.kind, sticks=(e, w)))
                elif d.mode == "bins":
                    fw = max(opt.fwhm, 0.0)
                    grid = make_grid(e.min() - 6 * max(fw, opt.de), e.max() + 6 * max(fw, opt.de), opt.de)
                    curve = broaden(e, y, grid, fw, opt.shape, opt.mixing)
                    out.append(Series(d.label, d.spin, grid, curve, d.name, ch, d.kind, sticks=(e, y),
                                      note="原始 .dos 直方图（Density 为归一化态分数）"))
                else:
                    yy = y
                    if opt.occupied and "occupied" in d.channels and ch == "total":
                        yy = d.channels["occupied"]
                    if opt.rebroaden_curves and opt.fwhm > 0:
                        # extra Gaussian convolution of an already broadened curve
                        de = np.median(np.diff(e))
                        sig = opt.fwhm / 2.3548200450309493
                        half = int(math.ceil(5 * sig / de))
                        kx = np.arange(-half, half + 1) * de
                        ker = np.exp(-0.5 * (kx / sig) ** 2)
                        ker /= ker.sum()
                        yy = np.convolve(yy, ker, mode="same")
                    out.append(Series(d.label, d.spin, e, yy, d.name, ch, d.kind))
        if opt.merge and out:
            merged = {}
            for s in out:
                merged.setdefault((s.spin, s.channel), []).append(s)
            res = []
            for (sp, ch), lst in merged.items():
                x = _common_grid([s.x for s in lst])
                y = np.sum([np.interp(x, s.x, s.y, left=0.0, right=0.0) for s in lst], axis=0)
                sticks = None
                if all(s.sticks is not None for s in lst):
                    sticks = (np.concatenate([s.sticks[0] for s in lst]), np.concatenate([s.sticks[1] for s in lst]))
                name = "Σ(" + "+".join(dict.fromkeys(s.kind for s in lst)) + ")"
                res.append(Series(name + ("" if sp == "none" else f"[{sp}]"), sp, x, y, "merged", ch, "merged", sticks))
            out = res
        if opt.integrate:
            out = [Series(s.name, s.spin, s.x, cumulative_trapz(s.x, s.y), s.source, s.channel, s.kind, None,
                          "累计积分 (IDOS)") for s in out]
        return out, warn

    # ------------------------------------------------------------------ statistics
    @staticmethod
    def statistics(series: list[Series], emin=None, emax=None):
        rows = []
        for s in series:
            if s.sticks is not None:
                d = stick_descriptors(s.sticks[0], s.sticks[1], emin, emax)
                basis = "离散态"
            else:
                d = curve_descriptors(s.x, s.y, emin, emax)
                basis = "曲线"
            rows.append({"name": s.name, "channel": s.channel, "spin": s.spin, "basis": basis, **d.as_dict()})
        return rows


def _common_grid(xs):
    lo = min(x[0] for x in xs)
    hi = max(x[-1] for x in xs)
    de = min(float(np.median(np.diff(x))) for x in xs if x.size > 1)
    return make_grid(lo, hi, de)


def write_series_csv(path, series: list[Series], unit_factor=1.0):
    """Write all series on a common grid (interpolated, zero outside)."""
    if not series:
        raise ValueError("没有可导出的曲线")
    x = _common_grid([s.x for s in series])
    cols = [x / unit_factor]
    head = ["Energy"]
    for s in series:
        cols.append(np.interp(x, s.x, s.y, left=0.0, right=0.0) * unit_factor)
        head.append(f"{s.name}:{s.channel}")
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(head)
        for row in np.column_stack(cols):
            w.writerow([f"{v:.6g}" for v in row])


def write_stats_csv(path, rows):
    keys = ["name", "channel", "spin", "basis", "integral", "center", "width", "skewness", "kurtosis",
            "filling", "occ_integral", "occ_center", "upper_edge"]
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow(keys)
        for r in rows:
            w.writerow([r[k] if isinstance(r[k], str) else ("" if r[k] != r[k] else f"{r[k]:.6f}") for k in keys])
