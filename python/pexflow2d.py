# -*- coding: utf-8 -*-
"""
pexflow2d.py -- Polynomial Expansion Optical Flow Velocimetry (PExFlow 2D)
==========================================================================

Optimised, dependency-light implementation (NumPy + SciPy) of the method
described in

    M. Flores, F. Rinderknecht, T. Gallot, L. G. Sarasua, Y. Abraham,
    N. Barrere, D. Freire Caporale (2026). "Optical flow method based on
    polynomial expansion for particle-laden fluid velocimetry".
    Experiments in Fluids 67:72. https://doi.org/10.1007/s00348-026-04222-1

Algorithm
---------
1. Local polynomial expansion of arbitrary degree n (Okorn & Harguess 2014)
   of the image intensity, fitted by Gaussian-weighted least squares on a
   (2w+1)x(2w+1) window (separable correlations). Odd degrees are promoted
   to the next even degree (3 -> 4).
2. Generalised Farneback system, one equation per coefficient c_ab with
   1 <= a+b <= n-1:
        A  = (A1(x) + A2(x + d0)) / 2
        Db = -(c2(x + d0) - c1(x)) + A d0
        d  = (sum w A^T A)^-1 (sum w A^T Db),   w = Gaussian(sigma_int)
3. Coarse-to-fine pyramid (factor ``pyr_scale``) with the displacement of
   each level used as prior for the next one.
4. Normalised-convolution border treatment (no invented pixels at the
   image edges, particles that leave the field of view get zero weight).

Optimisations with respect to the v1.0 GUI core (identical results):
  * the polynomial expansion of every frame is computed ONCE per pyramid
    level and reused by the two image pairs that share it (~2x fewer
    expansions on videos);
  * border operators, pixel grids and confidence ramps are cached per
    image size, and the border correction is fully vectorised (einsum);
  * float32 arithmetic throughout (like OpenCV).

Output convention: ``dx`` = columns (+x to the right), ``dy`` = rows
(+y DOWNWARDS, image convention), in pixels. Use :func:`to_physical` to
obtain velocities with +V pointing upwards.

Quick start
-----------
>>> from pexflow2d import PExFlow2D
>>> est = PExFlow2D()                       # default parameters (v1.0 GUI)
>>> dx, dy = est(frame0, frame1)            # one image pair
>>> U, V = est.sequence(video, frame_gap=1) # (T-1, M, N) displacement fields

Command line
------------
    python pexflow2d.py video.mp4 --dt 0.01 --scale 1e-4 --out flow.npz
    python pexflow2d.py "frames/*.png" --window 3 --degree 2 --out flow.mat

Author: Maicol Flores (Instituto de Fisica, Facultad de Ciencias, UdelaR).
License: MIT (see LICENSE in the repository).
"""
from __future__ import annotations

import os
import time
import numpy as np
from scipy.ndimage import correlate1d as _sp_correlate1d, map_coordinates

try:                                    # optional, multi-threaded separable filters
    if os.environ.get("PEXFLOW_NO_OPENCV"):
        raise ImportError("disabled by PEXFLOW_NO_OPENCV")
    import cv2 as _cv2
    _CV_BORDER = {"constant": _cv2.BORDER_CONSTANT, "nearest": _cv2.BORDER_REPLICATE,
                  "mirror": _cv2.BORDER_REFLECT_101}
except Exception:                       # OpenCV not installed -> pure NumPy/SciPy
    _cv2 = None

__version__ = "1.1.0"
__all__ = ["PExFlow2D", "pexflow2d", "to_physical", "derived_fields"]

_ONE = {}


def correlate1d(img, k, axis, mode="nearest", cval=0.0):
    """1-D correlation along ``axis`` (SciPy semantics). Uses OpenCV when it
    is available (same result, several times faster on multi-core CPUs)."""
    if _cv2 is not None and img.dtype == np.float32 and img.ndim == 2 \
            and min(img.shape) > 1 and mode in _CV_BORDER:
        k = np.asarray(k, np.float32)
        one = _ONE.setdefault("f32", np.ones(1, np.float32))
        kx, ky = (k, one) if axis == 1 else (one, k)
        return _cv2.sepFilter2D(img, -1, kx, ky, borderType=_CV_BORDER[mode])
    return _sp_correlate1d(img, k, axis=axis, mode=mode, cval=cval)


# ---------------------------------------------------------------------------
# Polynomial-expansion operator
# ---------------------------------------------------------------------------
def _build_operator(w, deg, dtype):
    """Precompute everything that only depends on (window, degree)."""
    ya, xb = [], []
    for d in range(deg + 1):
        for a in range(d, -1, -1):
            ya.append(a)
            xb.append(d - a)
    ya = np.array(ya)
    xb = np.array(xb)
    nT = ya.size

    r = np.arange(-w, w + 1, dtype=float)
    u = r / w
    app = np.exp(-r ** 2 / (2 * (0.3 * w) ** 2))          # applicability

    # Moments of the 1-D applicability: m[p] = sum app * u^p over a range
    lut = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(ya, xb))}
    eqs, used = [], set()
    for i in range(nT):
        a, b = int(ya[i]), int(xb[i])
        if a + b < 1 or a + b > deg - 1:
            continue
        iy, ix = lut[(a + 1, b)], lut[(a, b + 1)]
        eqs.append((i, iy, a + 1, ix, b + 1))
        used.update((i, iy, ix))

    k = np.stack([app * u ** p for p in range(deg + 1)], axis=1)   # (2w+1, deg+1)
    op = dict(w=w, deg=deg, ya=ya, xb=xb, nT=nT, app=app, u=u, dtype=dtype,
              k=k.astype(dtype), eqs=eqs, used=np.array(sorted(used)),
              scale=(1.0 / float(w) ** (ya + xb)).astype(dtype))
    op["C"] = _class_operator(op, -w, w, -w, w)
    return op


def _class_operator(op, rlo, rhi, clo, chi):
    """Least-squares operator for a window truncated to rows [rlo, rhi] and
    columns [clo, chi] (normalised convolution, certainty 0 outside)."""
    w, app, u = op["w"], op["app"], op["u"]
    t = np.arange(-w, w + 1)
    p2 = np.arange(2 * op["deg"] + 1)
    wy = app * ((t >= rlo) & (t <= rhi))
    wx = app * ((t >= clo) & (t <= chi))
    my = (wy[None, :] * u[None, :] ** p2[:, None]).sum(1)
    mx = (wx[None, :] * u[None, :] ** p2[:, None]).sum(1)
    ya, xb = op["ya"], op["xb"]
    G = my[ya[:, None] + ya[None, :]] * mx[xb[:, None] + xb[None, :]]
    return np.linalg.pinv(G).astype(op["dtype"])


def _border_classes(L, w):
    idx = np.arange(L)
    return np.maximum(-w, -idx), np.minimum(w, L - 1 - idx)


class _ShapeCache:
    """Per-image-size quantities (border operators, grids, ramps)."""

    def __init__(self, op, H, W):
        w, dt = op["w"], op["dtype"]
        self.H, self.W = H, W
        rlo, rhi = _border_classes(H, w)
        clo, chi = _border_classes(W, w)
        self.rb = np.where((rlo > -w) | (rhi < w))[0]          # border rows
        self.cb = np.where((clo > -w) | (chi < w))[0]          # border columns
        self.rows_in = np.setdiff1d(np.arange(H), self.rb)
        # Operators for border rows: one per (row, column-class)
        ccls = {}
        cidx = np.empty(W, np.intp)
        for c in range(W):
            cidx[c] = ccls.setdefault((clo[c], chi[c]), len(ccls))
        ckeys = list(ccls.keys())
        self.row_ops = np.stack([
            np.stack([_class_operator(op, rlo[r], rhi[r], cl, ch) for (cl, ch) in ckeys])
            for r in self.rb]) if self.rb.size else None        # (nrb, ncls, nT, nT)
        self.row_cidx = cidx
        self.col_ops = np.stack([_class_operator(op, -w, w, clo[c], chi[c])
                                 for c in self.cb]) if self.cb.size else None
        Y, X = np.mgrid[0:H, 0:W]
        self.Y = Y.astype(dt)
        self.X = X.astype(dt)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _gauss_kernel(sigma):
    r = max(1, int(np.ceil(3 * sigma)))
    x = np.arange(-r, r + 1, dtype=float)
    k = np.exp(-x ** 2 / (2 * sigma ** 2))
    return k / k.sum()


def _smooth2(img, k, mode="nearest"):
    if _cv2 is not None and img.dtype == np.float32 and min(img.shape) > 1:
        k = np.asarray(k, np.float32)
        return _cv2.sepFilter2D(img, -1, k, k, borderType=_CV_BORDER[mode])
    out = _sp_correlate1d(img, k, axis=0, mode=mode, cval=0.0)
    return _sp_correlate1d(out, k, axis=1, mode=mode, cval=0.0)


def _pyramid_sizes(M, N, L, min_size, scale):
    sizes = [(M, N)]
    s = 1.0
    for _ in range(1, L):
        s *= scale
        sz = (int(round(M * s)), int(round(N * s)))
        if min(sz) < min_size or sz == sizes[0]:
            break
        sizes.insert(0, sz)
    return sizes


def _area_matrix(n_in, n_out):
    sc = n_in / n_out
    Wm = np.zeros((n_out, n_in))
    for d in range(n_out):
        a, b = d * sc, (d + 1) * sc
        for i in range(int(np.floor(a)), min(int(np.ceil(b)), n_in)):
            Wm[d, i] = min(b, i + 1) - max(a, i)
        Wm[d] /= Wm[d].sum()
    return Wm


def _resize_area(img, sz):
    if img.shape == tuple(sz):
        return img
    Wy = _area_matrix(img.shape[0], sz[0]).astype(img.dtype)
    Wx = _area_matrix(img.shape[1], sz[1]).astype(img.dtype)
    return Wy @ img @ Wx.T


def _resize_bilinear(img, sz):
    H, W = img.shape
    if (H, W) == tuple(sz):
        return img
    if _cv2 is not None and img.dtype == np.float32:
        return _cv2.resize(img, (sz[1], sz[0]), interpolation=_cv2.INTER_LINEAR)
    yc = np.clip((np.arange(sz[0]) + 0.5) * H / sz[0] - 0.5, 0, H - 1)
    xc = np.clip((np.arange(sz[1]) + 0.5) * W / sz[1] - 0.5, 0, W - 1)
    Yq, Xq = np.meshgrid(yc, xc, indexing="ij")
    return map_coordinates(img, [Yq, Xq], order=1, mode="nearest").astype(img.dtype)


def _gaussian_presmooth(img, sigma):
    if sigma <= 0:
        return img
    ks = max(int(np.round(sigma * 5)) | 1, 3)
    r = ks // 2
    x = np.arange(-r, r + 1, dtype=float)
    k = np.exp(-x ** 2 / (2 * sigma ** 2))
    k = (k / k.sum()).astype(img.dtype)
    out = correlate1d(img, k, axis=1, mode="mirror")
    return correlate1d(out, k, axis=0, mode="mirror")


def _fine_blur(img):
    k = np.array([0.25, 0.5, 0.25], dtype=img.dtype)
    return correlate1d(correlate1d(img, k, axis=1, mode="mirror"), k, axis=0, mode="mirror")


def _build_pyramid(img, sizes, presmooth):
    N = img.shape[1]
    nL = len(sizes)
    out = [None] * nL
    if presmooth:
        for l, sz in enumerate(sizes):
            s = sz[1] / N
            out[l] = img if l == nL - 1 else _resize_bilinear(
                _gaussian_presmooth(img, (1.0 / s - 1.0) * 0.5), sz)
    else:
        out[-1] = img
        for l in range(nL - 2, -1, -1):
            out[l] = _resize_area(out[l + 1], sizes[l])
    return out


def _upsample_flow(dx0, dy0, sz):
    h0, w0 = dx0.shape
    sy, sx = sz[0] / h0, sz[1] / w0
    if _cv2 is not None and dx0.dtype == np.float32:
        return (_cv2.resize(dx0, (sz[1], sz[0]), interpolation=_cv2.INTER_LINEAR) * np.float32(sx),
                _cv2.resize(dy0, (sz[1], sz[0]), interpolation=_cv2.INTER_LINEAR) * np.float32(sy))
    yc = np.clip((np.arange(sz[0]) + 0.5) / sy - 0.5, 0, h0 - 1)
    xc = np.clip((np.arange(sz[1]) + 0.5) / sx - 0.5, 0, w0 - 1)
    Yq, Xq = np.meshgrid(yc, xc, indexing="ij")
    return (map_coordinates(dx0, [Yq, Xq], order=1, mode="nearest") * sx,
            map_coordinates(dy0, [Yq, Xq], order=1, mode="nearest") * sy)


# ---------------------------------------------------------------------------
# Core numerical kernels
# ---------------------------------------------------------------------------
def _poly_expand(img, op, sc, border):
    """Local polynomial coefficients as contiguous planes (nT, H, W), in
    pixel units."""
    H, W = img.shape
    deg, nT = op["deg"], op["nT"]
    mode = "constant" if border == "normalized" else "nearest"
    Vy = [correlate1d(img, op["k"][:, a], axis=0, mode=mode) for a in range(deg + 1)]
    P = np.empty((nT, H, W), img.dtype)
    for i in range(nT):
        P[i] = correlate1d(Vy[op["ya"][i]], op["k"][:, op["xb"][i]], axis=1, mode=mode)
    coef = (op["C"] @ P.reshape(nT, -1)).reshape(nT, H, W)
    if border == "normalized":
        if sc.rb.size:   # border rows (all columns): operator depends on (row, col-class)
            Cs = sc.row_ops[:, sc.row_cidx]                      # (nrb, W, nT, nT)
            coef[:, sc.rb] = np.einsum("rcij,jrc->irc", Cs, P[:, sc.rb], optimize=True)
        if sc.cb.size:   # border columns (interior rows)
            Pc = P[:, sc.rows_in][:, :, sc.cb]                   # (nT, nri, ncb)
            coef[:, sc.rows_in[:, None], sc.cb[None, :]] = np.einsum(
                "cij,jrc->irc", sc.col_ops, Pc, optimize=True)
    coef *= op["scale"][:, None, None]
    return coef


def _solve_level(c1, c2, dxp, dyp, op, sc, wk, reg, border):
    H, W = dxp.shape
    dt = dxp.dtype
    yb = sc.Y + dyp
    xb = sc.X + dxp

    # bilinear warp of the coefficient planes of the second image: the flat
    # indices and weights are computed once and reused for every plane
    yc = np.clip(yb, 0, H - 1)
    xc = np.clip(xb, 0, W - 1)
    y0 = np.minimum(yc.astype(np.intp), max(H - 2, 0))
    x0 = np.minimum(xc.astype(np.intp), max(W - 2, 0))
    wy = (yc - y0).astype(dt)
    wx = (xc - x0).astype(dt)
    i00 = (y0 * W + x0).ravel()
    i01 = i00 + (W > 1)
    i10 = i00 + W * (H > 1)
    i11 = i10 + (W > 1)
    w11 = wy * wx
    w10 = wy - w11
    w01 = wx - w11
    w00 = 1 - wy - w01
    c2w = [None] * c2.shape[0]
    for p in op["used"]:
        q = c2[p].reshape(-1)
        c2w[p] = (w00 * q[i00].reshape(H, W) + w01 * q[i01].reshape(H, W) +
                  w10 * q[i10].reshape(H, W) + w11 * q[i11].reshape(H, W))

    if border == "normalized":
        valid = ((xb >= 0) & (xb <= W - 1) & (yb >= 0) & (yb <= H - 1)).astype(dt)
        R = op["w"] + 1.0
        dxe = np.minimum(np.minimum(sc.X, xb), np.minimum(W - 1 - sc.X, W - 1 - xb))
        dye = np.minimum(np.minimum(sc.Y, yb), np.minimum(H - 1 - sc.Y, H - 1 - yb))
        rx = np.clip((dxe + 1) / R, 0, 1)
        ry = np.clip((dye + 1) / R, 0, 1)
        cw = valid * (rx * ry) ** 4
        smode = "constant"
    else:
        cw = None
        smode = "nearest"

    G11 = np.zeros((H, W), dt); G12 = np.zeros((H, W), dt); G22 = np.zeros((H, W), dt)
    h1 = np.zeros((H, W), dt); h2 = np.zeros((H, W), dt)
    for i0, iy, fy, ix, fx in op["eqs"]:
        Ay = (0.5 * fy) * (c1[iy] + c2w[iy])
        Ax = (0.5 * fx) * (c1[ix] + c2w[ix])
        Db = c1[i0] - c2w[i0] + Ay * dyp + Ax * dxp
        G11 += Ay * Ay; G12 += Ay * Ax; G22 += Ax * Ax
        h1 += Ay * Db;  h2 += Ax * Db

    if cw is not None:
        G11 *= cw; G12 *= cw; G22 *= cw; h1 *= cw; h2 *= cw
    G11 = _smooth2(G11, wk, smode); G12 = _smooth2(G12, wk, smode)
    G22 = _smooth2(G22, wk, smode)
    h1 = _smooth2(h1, wk, smode);   h2 = _smooth2(h2, wk, smode)

    lam = reg * np.mean(G11 + G22) + np.finfo(dt).tiny
    G11 += lam; G22 += lam
    det = G11 * G22 - G12 ** 2
    dy = (G22 * h1 - G12 * h2) / det
    dx = (G11 * h2 - G12 * h1) / det

    if border == "normalized":
        alpha = np.clip(_smooth2(valid, wk, "nearest"), 0.0, 1.0)
        if alpha.min() < 1.0:                       # only when something is invalid
            kb = _gauss_kernel(4.0).astype(dt)
            den = _smooth2(alpha, kb, "constant") + 1e-12
            fxv = _smooth2(alpha * dx, kb, "constant") / den
            fyv = _smooth2(alpha * dy, kb, "constant") / den
            dx = alpha * dx + (1 - alpha) * fxv
            dy = alpha * dy + (1 - alpha) * fyv
    return dx, dy


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
class PExFlow2D:
    """PExFlow 2D estimator.

    Parameters
    ----------
    window_size : int   Radius w of the polynomial window (window = 2w+1 px).
    poly_degree : int   Degree of the local polynomial (>= 2; odd -> next even).
    pyr_levels  : int   Maximum number of pyramid levels (coarse-to-fine).
    iterations  : int   Refinements per pyramid level.
    sigma_int   : float Std [px] of the Gaussian neighbourhood integration.
    reg         : float Tikhonov regularisation (relative).
    border      : 'normalized' (recommended) or 'replicate'.
    pyr_scale   : float Scale factor between pyramid levels, in (0, 1).
    presmooth   : bool  OpenCV-style Gaussian pre-smoothing of each level.
    min_size    : int   Minimum side [px] of the coarsest level (None = auto).
    fine_blur   : bool  [0.25 0.5 0.25] blur of the finest level.
    precision   : 'float32' (fast) or 'float64'.
    """

    def __init__(self, window_size=3, poly_degree=2, pyr_levels=6, iterations=3,
                 sigma_int=2.0, reg=1e-6, border="normalized", pyr_scale=0.65,
                 presmooth=True, min_size=None, fine_blur=False, precision="float32"):
        if int(poly_degree) != poly_degree or poly_degree < 2:
            raise ValueError("poly_degree must be an integer >= 2.")
        if int(window_size) != window_size or window_size < 1:
            raise ValueError("window_size must be an integer >= 1.")
        if not 0 < pyr_scale < 1:
            raise ValueError("pyr_scale must be in (0, 1).")
        if border not in ("normalized", "replicate"):
            raise ValueError("border must be 'normalized' or 'replicate'.")
        if precision not in ("float32", "float64"):
            raise ValueError("precision must be 'float32' or 'float64'.")
        poly_degree = int(poly_degree) + int(poly_degree) % 2      # parity correction
        self.params = dict(window_size=int(window_size), poly_degree=poly_degree,
                           pyr_levels=int(pyr_levels), iterations=int(iterations),
                           sigma_int=float(sigma_int), reg=float(reg), border=border,
                           pyr_scale=float(pyr_scale), presmooth=bool(presmooth),
                           min_size=min_size, fine_blur=bool(fine_blur), precision=precision)
        self.w = int(window_size)
        self.dtype = np.dtype(precision)
        self.op = _build_operator(self.w, poly_degree, self.dtype)
        self.wk = _gauss_kernel(sigma_int).astype(self.dtype)
        self._shape_cache = {}

    def __repr__(self):
        return "PExFlow2D(" + ", ".join(f"{k}={v!r}" for k, v in self.params.items()) + ")"

    # -- internals ---------------------------------------------------------
    def _sizes(self, M, N):
        p = self.params
        ms = max(2 * self.w + 1, 8) if p["min_size"] is None else max(int(p["min_size"]), 2 * self.w + 1)
        return _pyramid_sizes(M, N, p["pyr_levels"], ms, p["pyr_scale"])

    def _sc(self, shape):
        sc = self._shape_cache.get(shape)
        if sc is None:
            sc = self._shape_cache[shape] = _ShapeCache(self.op, *shape)
        return sc

    def _expand(self, img, sizes):
        """Pyramid + polynomial expansion of one frame (list over levels)."""
        f = np.asarray(img, dtype=self.dtype)
        if self.params["fine_blur"]:
            f = _fine_blur(f)
        pyr = _build_pyramid(f, sizes, self.params["presmooth"])
        b = self.params["border"]
        return [_poly_expand(I, self.op, self._sc(I.shape), b) for I in pyr]

    def _solve(self, E1, E2, sizes):
        p = self.params
        dx = np.zeros(sizes[0], self.dtype)
        dy = np.zeros(sizes[0], self.dtype)
        for l, sz in enumerate(sizes):
            if l > 0:
                dx, dy = _upsample_flow(dx, dy, sz)
            sc = self._sc(sz)
            for _ in range(max(1, p["iterations"])):
                dx, dy = _solve_level(E1[l], E2[l], dx, dy, self.op, sc, self.wk,
                                      p["reg"], p["border"])
        return dx, dy

    # -- public ------------------------------------------------------------
    def __call__(self, img1, img2):
        """Displacement (dx, dy) [px] from ``img1`` to ``img2`` (2-D arrays)."""
        img1 = np.asarray(img1)
        if img1.ndim != 2 or np.shape(img2) != img1.shape:
            raise ValueError("img1 and img2 must be 2-D arrays of the same shape.")
        sizes = self._sizes(*img1.shape)
        return self._solve(self._expand(img1, sizes), self._expand(img2, sizes), sizes)

    def sequence(self, frames, frame_gap=1, progress=False):
        """Displacement fields for a whole sequence.

        frames : array (T, M, N) or iterable of 2-D frames.
        Returns dX, dY with shape (T - frame_gap, M, N) [px]. Each frame is
        expanded once and reused by every pair that contains it.
        """
        frames = [np.asarray(f) for f in frames]
        T = len(frames)
        if T <= frame_gap:
            raise ValueError("The sequence must have more frames than frame_gap.")
        M, N = frames[0].shape
        sizes = self._sizes(M, N)
        n = T - frame_gap
        dX = np.empty((n, M, N), self.dtype)
        dY = np.empty((n, M, N), self.dtype)
        cache = {}
        t0 = time.perf_counter()
        for t in range(n):
            for k in (t, t + frame_gap):
                if k not in cache:
                    cache[k] = self._expand(frames[k], sizes)
            dX[t], dY[t] = self._solve(cache[t], cache[t + frame_gap], sizes)
            cache.pop(t, None)                      # frame t is no longer needed
            if progress:
                el = time.perf_counter() - t0
                print(f"\r  PExFlow2D  pair {t + 1}/{n}  ({el:.1f} s, "
                      f"{el / (t + 1):.2f} s/pair)", end="", flush=True)
        if progress:
            print()
        return dX, dY


def pexflow2d(video, window_size=3, poly_degree=2, pyr_levels=6, iterations=3,
              sigma_int=2.0, frame_gap=1, reg=1e-6, border="normalized", **kw):
    """Same signature as the MATLAB/Octave ``PExFlow2D`` function.

    video : array (M, N, T). Returns dX, dY with shape (M, N, T - frame_gap).
    """
    video = np.asarray(video)
    est = PExFlow2D(window_size, poly_degree, pyr_levels, iterations, sigma_int,
                    reg, border, **kw)
    dX, dY = est.sequence(np.moveaxis(video, 2, 0), frame_gap)
    return np.moveaxis(dX, 0, 2), np.moveaxis(dY, 0, 2)


def to_physical(dx, dy, dt=1.0, scale=1.0):
    """Pixel displacements -> velocities [scale units / s] with +V upwards."""
    return dx * (scale / dt), -dy * (scale / dt)


def derived_fields(U, V, h=1.0):
    """Derived quantities of a physical velocity field (U right, V up).

    U, V : 2-D arrays, row 0 at the TOP of the image. h : grid spacing.
    Returns dict with speed, vorticity (dV/dx - dU/dy), divergence and angle.
    """
    dUdy_img, dUdx = np.gradient(U, h)
    dVdy_img, dVdx = np.gradient(V, h)
    dUdy, dVdy = -dUdy_img, -dVdy_img          # image rows grow downwards
    return dict(speed=np.hypot(U, V), vorticity=dVdx - dUdy,
                divergence=dUdx + dVdy, angle=np.degrees(np.arctan2(V, U)))


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------
def _load_frames(src):
    import glob
    import os
    files = sorted(glob.glob(src))
    if len(files) > 1 or (files and os.path.splitext(files[0])[1].lower()
                          in (".png", ".tif", ".tiff", ".jpg", ".jpeg", ".bmp")):
        import imageio.v3 as iio
        frames = [iio.imread(f) for f in files]
    else:
        import cv2                                   # video file
        cap = cv2.VideoCapture(src)
        frames = []
        while True:
            ok, fr = cap.read()
            if not ok:
                break
            frames.append(fr)
        cap.release()
    if not frames:
        raise FileNotFoundError(f"No frames found in {src!r}")
    out = []
    for f in frames:
        f = np.asarray(f, np.float32)
        out.append(f.mean(axis=2) if f.ndim == 3 else f)
    return out


def _main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="PExFlow 2D - polynomial expansion optical flow velocimetry")
    ap.add_argument("source", help="video file or glob of images, e.g. 'frames/*.png'")
    ap.add_argument("--out", default="pexflow_result.npz", help=".npz or .mat output")
    ap.add_argument("--dt", type=float, default=1.0, help="time between frames [s]")
    ap.add_argument("--scale", type=float, default=1.0, help="calibration [m/px]")
    ap.add_argument("--gap", type=int, default=1, help="frame gap K")
    ap.add_argument("--window", type=int, default=3)
    ap.add_argument("--degree", type=int, default=2)
    ap.add_argument("--levels", type=int, default=6)
    ap.add_argument("--iterations", type=int, default=3)
    ap.add_argument("--sigma", type=float, default=2.0)
    a = ap.parse_args(argv)

    frames = _load_frames(a.source)
    est = PExFlow2D(a.window, a.degree, a.levels, a.iterations, a.sigma)
    print(f"{len(frames)} frames of {frames[0].shape[1]}x{frames[0].shape[0]} px  |  {est}")
    dX, dY = est.sequence(frames, a.gap, progress=True)
    U, V = to_physical(dX, dY, a.dt * a.gap, a.scale)
    M, N = frames[0].shape
    res = dict(U=np.moveaxis(U, 0, 2), V=np.moveaxis(V, 0, 2),
               x=np.arange(N) * a.scale, y=np.arange(M) * a.scale,
               t=np.arange(U.shape[0]) * a.dt,
               convention_V="positive = upwards (decreasing pixel row)")
    if a.out.lower().endswith(".mat"):
        from scipy.io import savemat
        savemat(a.out, res, do_compression=True)
    else:
        np.savez_compressed(a.out, **res)
    print(f"Saved {a.out}")


if __name__ == "__main__":
    _main()
