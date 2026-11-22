"""Matplotlib drawing helpers shared by the GUI and the command line (no Qt dependency)."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Validated categorical order (fixed, never cycled beyond 8 -> greyed "others")
PALETTE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
OVERFLOW = "#8c8c8c"
INK = "#222222"
MUTED = "#6b6b6b"
GRID = "#e6e6e6"

CJK_FONTS = ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "Source Han Sans SC", "WenQuanYi Micro Hei",
             "Arial Unicode MS", "DejaVu Sans"]


@dataclass
class PlotStyle:
    font_size: float = 11.0
    line_width: float = 1.6
    fill: bool = False
    fill_alpha: float = 0.25
    legend: bool = True
    grid: bool = True
    fermi_line: bool = True
    ticks_in: bool = True
    vertical: bool = False         # energy on the y axis (LOBSTER style)
    mirror_beta: bool = True
    xmin: float | None = None      # energy window
    xmax: float | None = None
    ymin: float | None = None
    ymax: float | None = None
    energy_label: str = "E − E$_F$ (eV)"
    title: str = ""
    shade_bonding: bool = True


def apply_rc():
    import matplotlib as mpl
    mpl.rcParams["font.sans-serif"] = CJK_FONTS + list(mpl.rcParams.get("font.sans-serif", []))
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["axes.unicode_minus"] = False
    mpl.rcParams["svg.fonttype"] = "none"
    mpl.rcParams["pdf.fonttype"] = 42


def color(i):
    return PALETTE[i] if i < len(PALETTE) else OVERFLOW


def _style_axes(ax, st: PlotStyle):
    ax.tick_params(direction="in" if st.ticks_in else "out", labelsize=st.font_size * 0.9, colors=INK,
                   top=st.ticks_in, right=st.ticks_in)
    for sp in ax.spines.values():
        sp.set_color("#555555")
        sp.set_linewidth(0.9)
    if st.grid:
        ax.grid(True, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)


def _autoscale(ax, data, st, vertical, sign=1.0, mirror=False):
    """Scale the value axis to the data inside the energy window."""
    lo, hi = st.xmin, st.xmax
    vmin, vmax = 0.0, 0.0
    for c in data:
        x = np.asarray(c["x"])
        y = sign * np.asarray(c["y"])
        if mirror and c.get("spin") == "beta":
            y = -y
        m = np.ones(x.size, bool)
        if lo is not None:
            m &= x >= lo
        if hi is not None:
            m &= x <= hi
        if m.any():
            vmin = min(vmin, float(y[m].min()))
            vmax = max(vmax, float(y[m].max()))
    if vmax == vmin:
        return
    pad = 0.06 * (vmax - vmin)
    lim = (vmin - pad if vmin < 0 else 0.0, vmax + pad if vmax > 0 else 0.0)
    (ax.set_xlim if vertical else ax.set_ylim)(*lim)


def _limits(ax, st, vertical):
    if vertical:
        if st.xmin is not None or st.xmax is not None:
            ax.set_ylim(st.xmin, st.xmax)
        if st.ymin is not None or st.ymax is not None:
            ax.set_xlim(st.ymin, st.ymax)
    else:
        if st.xmin is not None or st.xmax is not None:
            ax.set_xlim(st.xmin, st.xmax)
        if st.ymin is not None or st.ymax is not None:
            ax.set_ylim(st.ymin, st.ymax)


def draw_lines(fig, series, st: PlotStyle, ylabel="DOS (states/eV)", fermi=0.0):
    """series: list of dicts {label, x, y, spin, color_key}. Beta spin drawn negative if mirror_beta."""
    fig.clear()
    ax = fig.add_subplot(111)
    _style_axes(ax, st)
    keys = []
    for s in series:
        k = s.get("color_key", s["label"])
        if k not in keys:
            keys.append(k)
    for s in series:
        c = color(keys.index(s.get("color_key", s["label"])))
        y = np.asarray(s["y"])
        if st.mirror_beta and s.get("spin") == "beta":
            y = -y
        ls = "--" if (s.get("spin") == "beta" and not st.mirror_beta) else "-"
        lab = s["label"]
        if st.vertical:
            ax.plot(y, s["x"], color=c, lw=st.line_width, ls=ls, label=lab)
            if st.fill:
                ax.fill_betweenx(s["x"], 0, y, color=c, alpha=st.fill_alpha, lw=0)
        else:
            ax.plot(s["x"], y, color=c, lw=st.line_width, ls=ls, label=lab)
            if st.fill:
                ax.fill_between(s["x"], 0, y, color=c, alpha=st.fill_alpha, lw=0)
    if st.fermi_line and fermi is not None:
        (ax.axhline if st.vertical else ax.axvline)(fermi, color=MUTED, lw=0.9, ls="--")
    if any(s.get("spin") == "beta" for s in series) and st.mirror_beta:
        (ax.axvline if st.vertical else ax.axhline)(0, color="#999999", lw=0.7)
    if st.vertical:
        ax.set_ylabel(st.energy_label, fontsize=st.font_size, color=INK)
        ax.set_xlabel(ylabel, fontsize=st.font_size, color=INK)
    else:
        ax.set_xlabel(st.energy_label, fontsize=st.font_size, color=INK)
        ax.set_ylabel(ylabel, fontsize=st.font_size, color=INK)
    _autoscale(ax, series, st, st.vertical, mirror=st.mirror_beta)
    _limits(ax, st, st.vertical)
    if st.title:
        ax.set_title(st.title, fontsize=st.font_size * 1.05, color=INK)
    if st.legend and len({s['label'] for s in series}) >= 2:
        ax.legend(fontsize=st.font_size * 0.85, frameon=False, loc="best")
    fig.tight_layout()
    return ax


def draw_cohp(fig, curves, cumul, st: PlotStyle, quantity="COHP", negate=True, show_integral=True, fermi=0.0):
    """curves/cumul: list of dicts {label, x, y, spin}. LOBSTER-style: −COHP vs energy.

    The integrated curve goes into its own panel (never a second y scale on the same axes).
    """
    fig.clear()
    sign = -1.0 if (negate and quantity == "COHP") else 1.0
    qlab = ("−" if sign < 0 else "") + quantity
    n = 2 if (show_integral and cumul) else 1
    if st.vertical:
        axes = fig.subplots(1, n, sharey=True, squeeze=False)[0]
    else:
        axes = fig.subplots(n, 1, sharex=True, squeeze=False)[:, 0]
    labels = []
    for c in curves:
        if c["label"] not in labels:
            labels.append(c["label"])
    for pi, (ax, data, name) in enumerate(zip(axes, [curves, cumul][:n], [qlab, ("−I" + quantity) if sign < 0 else ("I" + quantity)])):
        _style_axes(ax, st)
        for c in data:
            col = color(labels.index(c["label"])) if c["label"] in labels else OVERFLOW
            y = sign * np.asarray(c["y"])
            ls = "--" if c.get("spin") == "down" else "-"
            lab = c["label"] + ("" if c.get("spin") in (None, "total") else f" ({c['spin']})")
            if st.vertical:
                ax.plot(y, c["x"], color=col, lw=st.line_width, ls=ls, label=lab)
                if st.fill and pi == 0:
                    ax.fill_betweenx(c["x"], 0, y, color=col, alpha=st.fill_alpha, lw=0)
            else:
                ax.plot(c["x"], y, color=col, lw=st.line_width, ls=ls, label=lab)
                if st.fill and pi == 0:
                    ax.fill_between(c["x"], 0, y, color=col, alpha=st.fill_alpha, lw=0)
        zero = ax.axvline if st.vertical else ax.axhline
        zero(0, color="#999999", lw=0.7)
        if st.fermi_line and fermi is not None:
            (ax.axhline if st.vertical else ax.axvline)(fermi, color=MUTED, lw=0.9, ls="--")
        if pi == 0 and st.shade_bonding and quantity in ("COHP", "COOP", "COBI"):
            # annotate which side is bonding
            bonding_positive = (sign < 0) or quantity in ("COOP", "COBI")
            txt_b, txt_a = ("bonding", "antibonding") if bonding_positive else ("antibonding", "bonding")
            if st.vertical:
                ax.text(0.98, 0.98, txt_b + " →", transform=ax.transAxes, ha="right", va="top",
                        fontsize=st.font_size * 0.75, color=MUTED)
                ax.text(0.02, 0.98, "← " + txt_a, transform=ax.transAxes, ha="left", va="top",
                        fontsize=st.font_size * 0.75, color=MUTED)
            else:
                ax.text(0.01, 0.97, "↑ " + txt_b, transform=ax.transAxes, ha="left", va="top",
                        fontsize=st.font_size * 0.75, color=MUTED)
                ax.text(0.01, 0.03, "↓ " + txt_a, transform=ax.transAxes, ha="left", va="bottom",
                        fontsize=st.font_size * 0.75, color=MUTED)
        unit = {"COHP": " (eV/eV)", "COOP": " (1/eV)", "COBI": " (1/eV)"}.get(quantity, "")
        iunit = {"COHP": " (eV)", "COOP": " (e)", "COBI": ""}.get(quantity, "")
        lab = name + (unit if pi == 0 else iunit)
        if st.vertical:
            ax.set_xlabel(lab, fontsize=st.font_size, color=INK)
            if pi == 0:
                ax.set_ylabel(st.energy_label, fontsize=st.font_size, color=INK)
        else:
            ax.set_ylabel(lab, fontsize=st.font_size, color=INK)
            if pi == n - 1:
                ax.set_xlabel(st.energy_label, fontsize=st.font_size, color=INK)
        _autoscale(ax, data, st, st.vertical, sign=sign)
        if st.vertical:
            if st.xmin is not None or st.xmax is not None:
                ax.set_ylim(st.xmin, st.xmax)
        else:
            if st.xmin is not None or st.xmax is not None:
                ax.set_xlim(st.xmin, st.xmax)
        if pi == 0 and (st.ymin is not None or st.ymax is not None):
            (ax.set_xlim if st.vertical else ax.set_ylim)(st.ymin, st.ymax)
        if pi == 0 and st.legend and len(data) >= 2:
            ax.legend(fontsize=st.font_size * 0.8, frameon=False, loc="best")
    if st.title:
        fig.suptitle(st.title, fontsize=st.font_size * 1.05, color=INK)
    fig.tight_layout()
    return axes
