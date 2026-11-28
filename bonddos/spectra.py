"""Broadening, integration and band-descriptor utilities (pure NumPy)."""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .constants import FWHM_PER_SIGMA

# NumPy 2.0 renamed trapz -> trapezoid (trapz was removed in 2.4); support both.
_trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz")
trapz = _trapz


def _erf_np(x):
    """erf with |error| < 1.2e-7 (Numerical Recipes erfc Chebyshev fit); NumPy only."""
    x = np.asarray(x, float)
    z = np.abs(x)
    t = 1.0 / (1.0 + 0.5 * z)
    ans = t * np.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 + t * (
        -0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))))
    return np.where(x >= 0, 1.0 - ans, ans - 1.0)


def _hilbert(y):
    """Analytic signal via FFT (same as scipy.signal.hilbert)."""
    y = np.asarray(y, float)
    n = y.size
    f = np.fft.fft(y)
    h = np.zeros(n)
    if n % 2 == 0:
        h[0] = h[n // 2] = 1
        h[1:n // 2] = 2
    else:
        h[0] = 1
        h[1:(n + 1) // 2] = 2
    return np.fft.ifft(f * h)


def make_grid(emin: float, emax: float, de: float) -> np.ndarray:
    if not (emax > emin):
        raise ValueError("能量上限必须大于下限")
    if de <= 0:
        raise ValueError("能量步长必须为正")
    n = int(round((emax - emin) / de)) + 1
    if n > 2_000_000:
        raise ValueError("能量网格点数过多，请增大步长或缩小能窗")
    return emin + de * np.arange(n)


def kernel(x: np.ndarray, fwhm: float, shape: str = "gaussian", mixing: float = 0.5) -> np.ndarray:
    """Unit-area line shape evaluated at x (same definitions as CP2K qs_dos_utils)."""
    shape = shape.lower()
    sigma = fwhm / FWHM_PER_SIGMA
    gamma = 0.5 * fwhm
    g = np.exp(-0.5 * (x / sigma) ** 2) / (sigma * math.sqrt(2 * math.pi))
    if shape.startswith("g"):
        return g
    lor = gamma / (math.pi * (x * x + gamma * gamma))
    if shape.startswith("l"):
        return lor
    eta = min(1.0, max(0.0, mixing))
    return eta * lor + (1 - eta) * g


def kernel_cdf(x: np.ndarray, fwhm: float, shape: str = "gaussian", mixing: float = 0.5) -> np.ndarray:
    shape = shape.lower()
    sigma = fwhm / FWHM_PER_SIGMA
    gamma = 0.5 * fwhm
    g = 0.5 * (1 + _erf_np(x / (sigma * math.sqrt(2))))
    if shape.startswith("g"):
        return g
    lor = 0.5 + np.arctan(x / gamma) / math.pi
    if shape.startswith("l"):
        return lor
    eta = min(1.0, max(0.0, mixing))
    return eta * lor + (1 - eta) * g


def broaden(energies, weights, grid, fwhm, shape="gaussian", mixing=0.5, cumulative=False, chunk=4000):
    """Sum of line shapes: out[e, c] = sum_n weights[n, c] K(grid[e] - energies[n]).

    weights may be 1-D (one channel) or 2-D (nstates, nchannels).
    With cumulative=True the running integral (CDF) is returned instead.
    """
    energies = np.asarray(energies, dtype=float).ravel()
    w = np.asarray(weights, dtype=float)
    one = w.ndim == 1
    if one:
        w = w[:, None]
    grid = np.asarray(grid, dtype=float)
    out = np.zeros((grid.size, w.shape[1]))
    if energies.size == 0:
        return out[:, 0] if one else out
    if fwhm <= 0:
        # stick spectrum histogram on the grid (bin = grid spacing), area preserving
        if grid.size < 2:
            raise ValueError("网格点过少")
        de = grid[1] - grid[0]
        idx = np.rint((energies - grid[0]) / de).astype(int)
        ok = (idx >= 0) & (idx < grid.size)
        for c in range(w.shape[1]):
            if cumulative:
                h = np.zeros(grid.size)
                np.add.at(h, idx[ok], w[ok, c])
                out[:, c] = np.cumsum(h) + w[(idx < 0), c].sum()
            else:
                np.add.at(out[:, c], idx[ok], w[ok, c] / de)
        return out[:, 0] if one else out
    cut = (8.0 if shape.lower().startswith("g") else 50.0) * fwhm
    order = np.argsort(energies)
    energies = energies[order]
    w = w[order]
    if cumulative:
        # states far below the grid contribute their full weight
        below = energies < grid[0] - cut
        out += w[below].sum(axis=0)[None, :]
        sel = ~below
        e2, w2 = energies[sel], w[sel]
        for start in range(0, e2.size, chunk):
            e = e2[start:start + chunk]
            k = kernel_cdf(grid[:, None] - e[None, :], fwhm, shape, mixing)
            out += k @ w2[start:start + chunk]
        return out[:, 0] if one else out
    lo = np.searchsorted(energies, grid[0] - cut)
    hi = np.searchsorted(energies, grid[-1] + cut)
    e_all, w_all = energies[lo:hi], w[lo:hi]
    # chunk over grid as well to cap memory (grid x chunk)
    gstep = max(1, int(2.0e7 // max(1, min(chunk, e_all.size))))
    for start in range(0, e_all.size, chunk):
        e = e_all[start:start + chunk]
        ww = w_all[start:start + chunk]
        for g0 in range(0, grid.size, gstep):
            gg = grid[g0:g0 + gstep]
            d = gg[:, None] - e[None, :]
            k = kernel(d, fwhm, shape, mixing)
            k[np.abs(d) > cut] = 0.0
            out[g0:g0 + gstep] += k @ ww
    return out[:, 0] if one else out


def cumulative_trapz(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    out = np.zeros_like(y)
    if y.ndim == 1:
        out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
    else:
        out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x)[:, None], axis=0)
    return out


def _window(x, y, a, b):
    """Restrict a piecewise-linear curve to [a, b] with interpolated end points."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    a = x[0] if a is None else max(a, x[0])
    b = x[-1] if b is None else min(b, x[-1])
    if b <= a:
        return np.array([a, b]), np.zeros(2)
    inside = (x > a) & (x < b)
    xs = np.concatenate([[a], x[inside], [b]])
    ys = np.concatenate([[np.interp(a, x, y)], y[inside], [np.interp(b, x, y)]])
    return xs, ys


@dataclass
class BandDescriptors:
    integral: float = float("nan")
    center: float = float("nan")
    width: float = float("nan")          # sqrt of 2nd central moment
    skewness: float = float("nan")
    kurtosis: float = float("nan")
    filling: float = float("nan")        # fraction of the window integral below E=0
    occ_integral: float = float("nan")   # integral below E=0 inside the window
    occ_center: float = float("nan")     # centre of the part below E=0
    upper_edge: float = float("nan")     # energy of the maximum of the Hilbert transform

    def as_dict(self):
        return dict(self.__dict__)


def curve_descriptors(x, y, emin=None, emax=None, zero=0.0) -> BandDescriptors:
    """Moments of a (non-negative) DOS curve inside [emin, emax]; energies relative to `zero`."""
    xs, ys = _window(np.asarray(x) - zero, y, None if emin is None else emin - 0.0, None if emax is None else emax)
    out = BandDescriptors()
    m0 = _trapz(ys, xs)
    out.integral = float(m0)
    if abs(m0) > 1e-14:
        c = _trapz(xs * ys, xs) / m0
        out.center = float(c)
        mu2 = _trapz((xs - c) ** 2 * ys, xs) / m0
        if mu2 > 0:
            out.width = float(math.sqrt(mu2))
            out.skewness = float(_trapz((xs - c) ** 3 * ys, xs) / m0 / mu2 ** 1.5)
            out.kurtosis = float(_trapz((xs - c) ** 4 * ys, xs) / m0 / mu2 ** 2)
    xo, yo = _window(xs, ys, None, 0.0)
    mo = _trapz(yo, xo) if xo[-1] > xo[0] else 0.0
    out.occ_integral = float(mo)
    if abs(m0) > 1e-14:
        out.filling = float(mo / m0)
    if abs(mo) > 1e-14:
        out.occ_center = float(_trapz(xo * yo, xo) / mo)
    if xs.size > 8 and abs(m0) > 1e-14:
        n = max(512, xs.size)
        xu = np.linspace(xs[0], xs[-1], n)
        yu = np.interp(xu, xs, ys)
        h = np.imag(_hilbert(yu))
        out.upper_edge = float(xu[int(np.argmax(h))])
    return out


def stick_descriptors(e, w, emin=None, emax=None, zero=0.0) -> BandDescriptors:
    e = np.asarray(e, float) - zero
    w = np.asarray(w, float)
    sel = np.ones(e.size, bool)
    if emin is not None:
        sel &= e >= emin
    if emax is not None:
        sel &= e <= emax
    e, w = e[sel], w[sel]
    out = BandDescriptors()
    m0 = w.sum()
    out.integral = float(m0)
    if abs(m0) > 1e-14:
        c = (e * w).sum() / m0
        out.center = float(c)
        mu2 = ((e - c) ** 2 * w).sum() / m0
        if mu2 > 0:
            out.width = float(math.sqrt(mu2))
            out.skewness = float(((e - c) ** 3 * w).sum() / m0 / mu2 ** 1.5)
            out.kurtosis = float(((e - c) ** 4 * w).sum() / m0 / mu2 ** 2)
    occ = e <= 0
    mo = w[occ].sum()
    out.occ_integral = float(mo)
    if abs(m0) > 1e-14:
        out.filling = float(mo / m0)
    if abs(mo) > 1e-14:
        out.occ_center = float((e[occ] * w[occ]).sum() / mo)
    return out
