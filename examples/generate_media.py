# -*- coding: utf-8 -*-
"""
generate_media.py -- reproduces every animation and figure of the README.

    pip install numpy scipy matplotlib pillow imageio imageio-ffmpeg
    python examples/generate_media.py            # all
    python examples/generate_media.py hero gyre  # only some

Outputs go to docs/media/ (animated WebP for the README + MP4 for download).
All flows are synthetic with an analytical ground truth (see synthetic.py).
"""
import os
import sys
import json
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from scipy.ndimage import map_coordinates

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "python"))
sys.path.insert(0, HERE)

import pexflow2d as pf                                         # noqa: E402
from pexflow2d import PExFlow2D, to_physical, derived_fields   # noqa: E402
from synthetic import FLOWS, ParticleMovie, r2_score, cosine_similarity  # noqa: E402

OUT = os.path.join(ROOT, "docs", "media")
os.makedirs(OUT, exist_ok=True)

BG, PANEL, FG, MUTED, GRID = "#0d1117", "#161b22", "#e6edf3", "#8b949e", "#30363d"
ACCENT, ACCENT2 = "#58a6ff", "#3fb950"


def style():
    plt.rcParams.update({
        "figure.facecolor": BG, "axes.facecolor": PANEL, "savefig.facecolor": BG,
        "text.color": FG, "axes.labelcolor": FG, "axes.edgecolor": GRID,
        "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlecolor": FG,
        "font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 11,
        "axes.titleweight": "bold", "legend.frameon": False,
    })


def to_rgb(fig):
    fig.canvas.draw()
    a = np.asarray(fig.canvas.buffer_rgba())[:, :, :3]
    h, w = a.shape[:2]
    return a[: h - h % 2, : w - w % 2].copy()


def save_anim(frames, name, fps=12, quality=82):
    """Animated WebP (shown in the README) + MP4 (for download)."""
    from PIL import Image
    imgs = [Image.fromarray(f) for f in frames]
    webp = os.path.join(OUT, name + ".webp")
    imgs[0].save(webp, format="WEBP", save_all=True, append_images=imgs[1:],
                 duration=int(1000 / fps), loop=0, quality=quality, method=6)
    try:
        import imageio.v2 as iio
        mp4 = os.path.join(OUT, name + ".mp4")
        with iio.get_writer(mp4, fps=fps, codec="libx264", quality=8,
                            pixelformat="yuv420p", macro_block_size=1) as w:
            for f in frames:
                w.append_data(f)
    except Exception as e:                                 # pragma: no cover
        print("  (mp4 skipped:", e, ")")
    print(f"  -> {webp}  ({os.path.getsize(webp) / 1e6:.1f} MB)")


def label(ax, text, loc=(0.012, 0.965)):
    ax.text(*loc, text, transform=ax.transAxes, ha="left", va="top", fontsize=9.5,
            color=FG, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", fc=BG, ec="none", alpha=0.72))


def clean(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def vorticity_px(dx, dy, smooth=0.0):
    U, V = to_physical(dx, dy)
    w = derived_fields(U, V)["vorticity"]
    if smooth:
        from scipy.ndimage import gaussian_filter
        w = gaussian_filter(w, smooth)
    return w


SMOOTH = 1.5        # display-only Gaussian filter of the vorticity maps [px]


# ---------------------------------------------------------------------------
# 1. Hero: von Karman-type vortex street
# ---------------------------------------------------------------------------
def hero(n=48):
    print("[hero] vortex street")
    H, W = 192, 480
    vel = FLOWS["vortex_street"](H, W)
    mov = ParticleMovie(vel, (H, W), ppp=0.05, inflow="left", seed=2)
    frames, (tx, ty) = mov.run(n + 1)
    est = PExFlow2D()
    t0 = time.perf_counter()
    dX, dY = est.sequence(frames, progress=True)
    tpp = (time.perf_counter() - t0) / n
    r2 = [r2_score(dX[t], dY[t], tx[t], ty[t], 4) for t in range(n)]
    vort = np.array([vorticity_px(dX[t], dY[t], SMOOTH) for t in range(n)])
    vmax = np.percentile(np.abs(vort), 99.5)
    smax = np.percentile(np.hypot(dX, dY), 99.8)
    s = 16
    ys, xs = np.mgrid[s // 2:H:s, s // 2:W:s]
    out = []
    for t in range(n):
        fig, axs = plt.subplots(3, 1, figsize=(6.6, 7.35), dpi=86)
        fig.subplots_adjust(0.01, 0.055, 0.99, 0.995, hspace=0.035)
        axs[0].imshow(frames[t], cmap="gray", vmin=0, vmax=200)
        label(axs[0], "Input: particle images")
        sp = np.hypot(dX[t], dY[t])
        axs[1].imshow(sp, cmap="magma", vmin=0, vmax=smax)
        axs[1].quiver(xs, ys, dX[t][ys, xs], -dY[t][ys, xs], color="w", scale=smax * 26,
                      width=0.0028, headwidth=3.6, alpha=0.9)
        label(axs[1], "PExFlow  |u|  + vectors")
        axs[2].imshow(vort[t], cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        label(axs[2], f"PExFlow vorticity  $\\omega_z$  (Gaussian $\\sigma$ = {SMOOTH} px)")
        for a in axs:
            clean(a)
        fig.text(0.5, 0.018, f"frame {t + 1:02d}/{n}    R$^2$ vs. ground truth = {r2[t]:.3f}"
                 f"    {tpp:.2f} s/frame ({W}x{H} px, CPU)", ha="center", color=MUTED, fontsize=9)
        out.append(to_rgb(fig))
        plt.close(fig)
    save_anim(out, "hero_vortex_street", fps=12)
    return dict(r2_mean=float(np.mean(r2)), s_per_pair=tpp)


# ---------------------------------------------------------------------------
# 2. Double gyre: streamlines + error against the analytical field
# ---------------------------------------------------------------------------
def gyre(n=48):
    print("[gyre] time-periodic double gyre")
    H, W = 200, 400
    vel = FLOWS["double_gyre"](H, W, period=n)
    mov = ParticleMovie(vel, (H, W), ppp=0.055, seed=4)
    frames, (tx, ty) = mov.run(n + 1)
    dX, dY = PExFlow2D().sequence(frames, progress=True)
    err = np.hypot(dX - tx, dY - ty)
    r2 = [r2_score(dX[t], dY[t], tx[t], ty[t], 4) for t in range(n)]
    smax = np.percentile(np.hypot(dX, dY), 99.5)
    emax = 0.5
    x, y = np.arange(W), np.arange(H)
    out = []
    for t in range(n):
        fig, axs = plt.subplots(3, 1, figsize=(5.6, 8.2), dpi=86)
        fig.subplots_adjust(0.01, 0.05, 0.99, 0.995, hspace=0.035)
        axs[0].imshow(frames[t], cmap="gray", vmin=0, vmax=200)
        label(axs[0], "Input: particle images")
        sp = np.hypot(dX[t], dY[t])
        axs[1].imshow(np.zeros((H, W)), cmap="gray", vmin=0, vmax=1)
        axs[1].streamplot(x, y, dX[t], dY[t], color=sp, cmap="viridis", norm=plt.Normalize(0, smax),
                          density=1.35, linewidth=0.9, arrowsize=0.7)
        axs[1].set_xlim(0, W - 1)
        axs[1].set_ylim(H - 1, 0)
        label(axs[1], "PExFlow streamlines (colour = |u|)")
        axs[2].imshow(err[t], cmap="inferno", vmin=0, vmax=emax)
        label(axs[2], f"|error| vs. analytical  (0 - {emax} px)")
        for a in axs:
            clean(a)
        fig.text(0.5, 0.016, f"t = {t + 1:02d}/{n}    R$^2$ = {r2[t]:.4f}    "
                 f"RMS error = {np.sqrt(np.mean(err[t] ** 2)):.3f} px",
                 ha="center", color=MUTED, fontsize=9)
        out.append(to_rgb(fig))
        plt.close(fig)
    save_anim(out, "double_gyre", fps=12)
    return dict(r2_mean=float(np.mean(r2)))


# ---------------------------------------------------------------------------
# 3. Lorentz-forced vortex lattice + Lagrangian tracers advected by PExFlow
# ---------------------------------------------------------------------------
def lorentz(n=48):
    print("[lorentz] vortex lattice + virtual tracers")
    H = W = 256
    vel = FLOWS["lorentz_cells"](H, W, U0=4.0, period=n)
    mov = ParticleMovie(vel, (H, W), ppp=0.06, periodic=True, seed=5)
    frames, _ = mov.run(n + 1, truth=False)
    dX, dY = PExFlow2D().sequence(frames, progress=True)
    vort = np.array([vorticity_px(dX[t], dY[t], SMOOTH) for t in range(n)])
    vmax = np.percentile(np.abs(vort), 99.5)
    s = 16
    ys, xs = np.mgrid[s // 2:H:s, s // 2:W:s]
    # virtual tracers seeded on a lattice, advected with the measured field
    g = np.linspace(14, W - 14, 13)
    px, py = [a.ravel() for a in np.meshgrid(g, g)]
    col = plt.cm.plasma((px + py) / (2 * W))
    trail = []
    out = []
    for t in range(n):
        trail.append((px.copy(), py.copy()))
        trail = trail[-14:]
        fig, axs = plt.subplots(1, 3, figsize=(11.2, 4.05), dpi=80)
        fig.subplots_adjust(0.005, 0.075, 0.995, 0.995, wspace=0.02)
        axs[0].imshow(frames[t], cmap="gray", vmin=0, vmax=200)
        label(axs[0], "Input: particle images")
        axs[1].imshow(vort[t], cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        axs[1].quiver(xs, ys, dX[t][ys, xs], -dY[t][ys, xs], color="k", scale=60,
                      width=0.004, alpha=0.75)
        label(axs[1], "PExFlow vorticity + vectors")
        axs[2].imshow(np.hypot(dX[t], dY[t]), cmap="Greys_r", vmin=0, vmax=12, alpha=1)
        for k, (qx, qy) in enumerate(trail):
            a = (k + 1) / len(trail)
            axs[2].scatter(qx, qy, s=4 + 10 * (k == len(trail) - 1), c=col, alpha=a * 0.9,
                           linewidths=0)
        axs[2].set_xlim(0, W - 1)
        axs[2].set_ylim(H - 1, 0)
        label(axs[2], "Virtual tracers advected by PExFlow")
        for a in axs:
            clean(a)
        fig.text(0.5, 0.022, f"Chessboard vortex lattice (Lorentz-forcing-like)    frame {t + 1:02d}/{n}",
                 ha="center", color=MUTED, fontsize=9)
        out.append(to_rgb(fig))
        plt.close(fig)
        # advect tracers (midpoint rule with the measured displacement)
        c = [np.clip(py, 0, H - 1), np.clip(px, 0, W - 1)]
        mx = map_coordinates(dX[t], c, order=1)
        my = map_coordinates(dY[t], c, order=1)
        c = [np.clip(py + my / 2, 0, H - 1), np.clip(px + mx / 2, 0, W - 1)]
        px = np.mod(px + map_coordinates(dX[t], c, order=1), W)
        py = np.mod(py + map_coordinates(dY[t], c, order=1), H)
    save_anim(out, "lorentz_lattice", fps=12)


# ---------------------------------------------------------------------------
# 4. Coarse-to-fine pipeline figure
# ---------------------------------------------------------------------------
def pipeline():
    print("[pipeline] coarse-to-fine figure")
    H = W = 256
    mov = ParticleMovie(FLOWS["lamb_oseen"](H, W), (H, W), seed=3)
    fr, (tx, ty) = mov.run(2)
    est = PExFlow2D()
    sizes = est._sizes(H, W)
    E1 = est._expand(fr[0], sizes)
    E2 = est._expand(fr[1], sizes)
    pyr = pf._build_pyramid(fr[0].astype(np.float32), sizes, True)
    dx = np.zeros(sizes[0], np.float32)
    dy = np.zeros(sizes[0], np.float32)
    rec = []
    for l, sz in enumerate(sizes):
        if l > 0:
            dx, dy = pf._upsample_flow(dx, dy, sz)
        sc = est._sc(sz)
        for _ in range(3):
            dx, dy = pf._solve_level(E1[l], E2[l], dx, dy, est.op, sc, est.wk, 1e-6, "normalized")
        rec.append((dx.copy(), dy.copy()))
    nL = len(sizes)
    fig, axs = plt.subplots(2, nL, figsize=(2.15 * nL, 4.75), dpi=120)
    fig.subplots_adjust(0.01, 0.02, 0.99, 0.86, wspace=0.06, hspace=0.08)
    for l in range(nL):
        h, w = sizes[l]
        axs[0, l].imshow(pyr[l], cmap="gray", vmin=0, vmax=150, interpolation="nearest")
        axs[0, l].set_title(f"L = {l + 1}   ({w}x{h})", fontsize=9.5, color=FG if l < nL - 1 else ACCENT2)
        f = W / w
        sp = np.hypot(*rec[l]) * f
        axs[1, l].imshow(sp, cmap="magma", vmin=0, vmax=8.5, interpolation="nearest")
        st = max(1, int(round(h / 11)))
        yy, xx = np.mgrid[st // 2:h:st, st // 2:w:st]
        axs[1, l].quiver(xx, yy, rec[l][0][yy, xx], -rec[l][1][yy, xx], color="w", width=0.012,
                         scale=6.5 * 8.5 / f, headwidth=3)
        r2 = r2_score(*[np.asarray(pf._upsample_flow(*rec[l], (H, W)))[i] for i in (0, 1)], tx[0], ty[0], 4)
        axs[1, l].text(0.5, -0.02, f"R$^2$ = {r2:.3f}", transform=axs[1, l].transAxes, ha="center",
                       va="top", color=MUTED, fontsize=9)
        for a in axs[:, l]:
            clean(a)
    fig.subplots_adjust(bottom=0.07)
    fig.text(0.5, 0.965, "Coarse-to-fine estimation: the displacement of each level is the prior of the next",
             ha="center", color=FG, fontsize=11.5, fontweight="bold")
    fig.text(0.01, 0.91, "coarse", color=MUTED, fontsize=9)
    fig.text(0.99, 0.91, "fine", color=MUTED, fontsize=9, ha="right")
    fig.add_artist(FancyArrowPatch((0.06, 0.915), (0.955, 0.915), transform=fig.transFigure,
                                   arrowstyle="-|>", mutation_scale=14, color=ACCENT, lw=1.4))
    fig.savefig(os.path.join(OUT, "pipeline.png"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# 5. Accuracy on canonical flows
# ---------------------------------------------------------------------------
CASES = [("Poiseuille", "poiseuille", {}), ("Lamb-Oseen vortex", "lamb_oseen", {}),
         ("Sink + swirl", "sink", {}), ("Vortex lattice", "lorentz_cells", {}),
         ("Double gyre", "double_gyre", {}), ("Vortex street", "vortex_street", {"inflow": "left"}),
         ("2-D turbulence", "turbulence", {"periodic": True})]


def accuracy(seeds=(1, 2, 3)):
    print("[accuracy] canonical flows")
    H = W = 256
    est = PExFlow2D()
    rows = []
    for title, key, kw in CASES:
        r2s, ss, rms, dmax = [], [], [], []
        for sd in seeds:
            mov = ParticleMovie(FLOWS[key](H, W), (H, W), seed=sd, **kw)
            for _ in range(5 if key == "vortex_street" else 0):
                mov.step(0)
            fr, (tx, ty) = mov.run(2)
            dx, dy = est(fr[0], fr[1])
            r2s.append(r2_score(dx, dy, tx[0], ty[0], 4))
            ss.append(cosine_similarity(dx, dy, tx[0], ty[0], 4))
            rms.append(float(np.sqrt(np.mean((dx - tx[0]) ** 2 + (dy - ty[0]) ** 2))))
            dmax.append(float(np.hypot(tx[0], ty[0]).max()))
        rows.append(dict(case=title, R2=float(np.mean(r2s)), S=float(np.mean(ss)),
                         rms_px=float(np.mean(rms)), max_disp_px=float(np.mean(dmax))))
        print(f"  {title:18s} R2={rows[-1]['R2']:.4f}  S={rows[-1]['S']:.4f}  RMS={rows[-1]['rms_px']:.3f}px")
    with open(os.path.join(OUT, "accuracy.json"), "w") as f:
        json.dump(rows, f, indent=1)

    fig, ax = plt.subplots(figsize=(9.2, 3.9), dpi=120)
    fig.subplots_adjust(0.07, 0.25, 0.99, 0.86)
    xk = np.arange(len(rows))
    b1 = ax.bar(xk - 0.19, [r["R2"] for r in rows], 0.36, color=ACCENT, label="R$^2$")
    b2 = ax.bar(xk + 0.19, [r["S"] for r in rows], 0.36, color=ACCENT2, label="cosine similarity S")
    for b in list(b1) + list(b2):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.004, f"{b.get_height():.3f}",
                ha="center", va="bottom", fontsize=7.5, color=FG, rotation=90)
    ax.set_xticks(xk)
    ax.set_xticklabels([r["case"] for r in rows], rotation=22, ha="right", color=FG)
    ax.set_ylim(0.8, 1.04)
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.legend(loc="upper right", ncol=2, bbox_to_anchor=(1.0, 1.19))
    ax.set_title("Accuracy vs. analytical ground truth  (defaults, 256x256 px)",
                 loc="left", fontsize=10.5, pad=10)
    fig.savefig(os.path.join(OUT, "accuracy.png"))
    plt.close(fig)
    return rows


# ---------------------------------------------------------------------------
# 6. Robustness to the parameters (window size x polynomial degree)
# ---------------------------------------------------------------------------
def robustness():
    print("[robustness] parameter map")
    H = W = 256
    flows = ["poiseuille", "lamb_oseen", "sink", "lorentz_cells", "turbulence"]
    data = {}
    for key in flows:
        mov = ParticleMovie(FLOWS[key](H, W), (H, W), seed=11, periodic=(key == "turbulence"))
        data[key] = mov.run(2)
    wins = list(range(1, 9))
    degs = [2, 4, 6, 8]
    R = np.zeros((len(degs), len(wins)))
    for i, d in enumerate(degs):
        for j, w in enumerate(wins):
            est = PExFlow2D(window_size=w, poly_degree=d)
            vals = []
            for key in flows:
                fr, (tx, ty) = data[key]
                dx, dy = est(fr[0], fr[1])
                vals.append(r2_score(dx, dy, tx[0], ty[0], 4))
            R[i, j] = np.mean(vals)
        print(f"  degree {d}: " + " ".join(f"{v:.3f}" for v in R[i]))
    np.save(os.path.join(OUT, "robustness.npy"), R)
    fig, ax = plt.subplots(figsize=(7.6, 3.5), dpi=120)
    fig.subplots_adjust(0.11, 0.17, 0.93, 0.85)
    im = ax.imshow(R, cmap="viridis", vmin=0.85, vmax=1.0, aspect="auto")
    for i in range(len(degs)):
        for j in range(len(wins)):
            ax.text(j, i, f"{R[i, j]:.3f}", ha="center", va="center", fontsize=8.5,
                    color="k" if R[i, j] > 0.95 else "w")
    ax.add_patch(plt.Rectangle((wins.index(3) - 0.5, degs.index(2) - 0.5), 1, 1, fill=False,
                               ec="#ff7b72", lw=2.2))
    ax.set_xticks(range(len(wins)))
    ax.set_xticklabels([f"{w}\n({2 * w + 1} px)" for w in wins], fontsize=8.5)
    ax.set_yticks(range(len(degs)))
    ax.set_yticklabels([f"{d - 1}-{d}" if d > 2 else "2" for d in degs])
    ax.set_xlabel("window_size  w   (window 2w+1)")
    ax.set_ylabel("poly_degree")
    cb = fig.colorbar(im, ax=ax, pad=0.015)
    cb.outline.set_visible(False)
    ax.set_title("Mean R$^2$ over 5 canonical flows -- red: default", fontsize=10.5)
    fig.savefig(os.path.join(OUT, "robustness.png"))
    plt.close(fig)


# ---------------------------------------------------------------------------
# 7. 2-D turbulence: vorticity + energy spectrum
# ---------------------------------------------------------------------------
def spectrum_1d(u, v):
    H, W = u.shape
    wy = np.hanning(H)[:, None] * np.hanning(W)[None, :]
    U = np.fft.fft2((u - u.mean()) * wy)
    V = np.fft.fft2((v - v.mean()) * wy)
    E2 = 0.5 * (np.abs(U) ** 2 + np.abs(V) ** 2)
    ky = np.fft.fftfreq(H) * H
    kx = np.fft.fftfreq(W) * W
    K = np.hypot(*np.meshgrid(kx, ky)).round().astype(int)
    E = np.bincount(K.ravel(), E2.ravel())
    return np.arange(E.size), E


def turbulence():
    print("[turbulence] vorticity + spectrum")
    H = W = 384
    vel = FLOWS["turbulence"](H, W, slope=-3.0, kmin=3, rms=1.5, seed=3)
    mov = ParticleMovie(vel, (H, W), ppp=0.07, periodic=True, seed=6)
    fr, (tx, ty) = mov.run(2)
    dx, dy = PExFlow2D()(fr[0], fr[1])
    from scipy.ndimage import gaussian_filter
    sg = 1.5     # same mild filter on both maps (display only)
    wt = gaussian_filter(vorticity_px(tx[0], ty[0]), sg)
    we = gaussian_filter(vorticity_px(dx, dy), sg)
    cc = np.corrcoef(wt[8:-8, 8:-8].ravel(), we[8:-8, 8:-8].ravel())[0, 1]
    vmax = np.percentile(np.abs(wt), 99.5)
    k, Et = spectrum_1d(tx[0], ty[0])
    _, Ee = spectrum_1d(dx, dy)
    fig = plt.figure(figsize=(12.4, 4.0), dpi=110)
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.25], left=0.01, right=0.985, bottom=0.13,
                          top=0.9, wspace=0.08)
    for i, (w, t) in enumerate([(wt, "Ground-truth vorticity"),
                                (we, f"PExFlow vorticity  (corr. {cc:.2f})")]):
        ax = fig.add_subplot(gs[i])
        ax.imshow(w, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        clean(ax)
        ax.set_title(t)
    fig.text(0.36, 0.06, f"maps filtered with the same Gaussian, $\\sigma$ = {sg} px", ha="center",
             color=MUTED, fontsize=9)
    ax = fig.add_subplot(gs[2])
    m = (k >= 1) & (k < H // 2)
    ax.loglog(k[m], Et[m], color=FG, lw=2, label="ground truth")
    ax.loglog(k[m], Ee[m], color=ACCENT, lw=1.8, ls="--", label="PExFlow")
    kr = np.array([6, 60])
    ax.loglog(kr, Et[6] * (kr / 6.0) ** -3 * 2.5, color="#ff7b72", lw=1.2, label="$k^{-3}$ (Kraichnan)")
    ax.set_xlabel("wavenumber k  [1/image]")
    ax.set_ylabel("E(k)")
    ax.grid(True, which="both", color=GRID, lw=0.4)
    ax.legend(loc="lower left")
    r2 = r2_score(dx, dy, tx[0], ty[0], 4)
    ax.set_title(f"Kinetic-energy spectrum   (R$^2$ = {r2:.3f})")
    fig.savefig(os.path.join(OUT, "turbulence_spectrum.png"))
    plt.close(fig)
    return dict(r2=float(r2), vort_corr=float(cc))


# ---------------------------------------------------------------------------
# 8. Banner
# ---------------------------------------------------------------------------
def banner():
    print("[banner]")
    H, W = 300, 1200
    vel = FLOWS["turbulence"](H, W, slope=-3.0, kmin=4, rms=1.4, seed=9)
    mov = ParticleMovie(vel, (H, W), ppp=0.05, periodic=True, seed=9)
    fr, _ = mov.run(2, truth=False)
    dx, dy = PExFlow2D()(fr[0], fr[1])
    w = vorticity_px(dx, dy)
    vmax = np.percentile(np.abs(w), 99)
    fig = plt.figure(figsize=(12, 3.0), dpi=110)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(w, cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.streamplot(np.arange(W), np.arange(H), dx, dy, color=(1, 1, 1, 0.45), density=2.0,
                  linewidth=0.6, arrowsize=0.5)
    ax.set_xlim(0, W - 1)
    ax.set_ylim(H - 1, 0)
    yy = np.linspace(-1, 1, H)[:, None] * np.ones((1, W))
    xx = np.ones((H, 1)) * np.linspace(-1, 1, W)[None, :]
    shade = np.clip(1.05 - 0.55 * (xx ** 2) - 0.25 * yy ** 2, 0, 1)
    rgba = np.zeros((H, W, 4))
    rgba[..., :3] = matplotlib.colors.to_rgb(BG)
    rgba[..., 3] = 0.25 + 0.5 * shade
    ax.imshow(rgba, extent=(-0.5, W - 0.5, H - 0.5, -0.5))
    clean(ax)
    ax.text(0.5, 0.60, "PExFlow", transform=ax.transAxes, ha="center", va="center",
            fontsize=58, fontweight="bold", color="white")
    ax.text(0.5, 0.33, "Polynomial Expansion Optical Flow Velocimetry", transform=ax.transAxes,
            ha="center", va="center", fontsize=19, color="#dfe6ee")
    ax.text(0.5, 0.17, "dense velocity fields from particle images  ·  Python · MATLAB · Octave",
            transform=ax.transAxes, ha="center", va="center", fontsize=11.5, color="#c9d1d9")
    fig.savefig(os.path.join(OUT, "banner.png"))
    plt.close(fig)


JOBS = dict(banner=banner, hero=hero, gyre=gyre, lorentz=lorentz, pipeline=pipeline,
            accuracy=accuracy, robustness=robustness, turbulence=turbulence)

if __name__ == "__main__":
    style()
    names = sys.argv[1:] or list(JOBS)
    summary = {}
    for nm in names:
        t0 = time.perf_counter()
        r = JOBS[nm]()
        if r is not None:
            summary[nm] = r
        print(f"  [{nm}] {time.perf_counter() - t0:.1f} s")
    if summary:
        p = os.path.join(OUT, "summary.json")
        old = json.load(open(p)) if os.path.exists(p) else {}
        old.update(summary)
        json.dump(old, open(p, "w"), indent=1)
