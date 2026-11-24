# -*- coding: utf-8 -*-
"""
benchmark.py -- execution time per image pair of the Python implementation.

    python examples/benchmark.py                       # this repository
    python examples/benchmark.py --v10 path/to/v1.0/pexflow2d.py --matlab 0.374

Measures the optimised core with and without OpenCV and, optionally, the
v1.0 core shipped with the GUI. MATLAB/Octave times can be added by hand
(run  matlab/benchmark_pexflow.m) and are plotted in docs/media/speed.png.
"""
import os
import sys
import json
import time
import argparse
import importlib.util

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "python"))
sys.path.insert(0, HERE)
import pexflow2d as pf                     # noqa: E402
from synthetic import FLOWS, ParticleMovie  # noqa: E402


def timeit(fn, reps=1):
    best = np.inf
    for _ in range(reps):
        t = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--size", type=int, default=512)
    ap.add_argument("--pairs", type=int, default=10)
    ap.add_argument("--v10", help="path to the v1.0 pexflow2d.py (GUI core) to compare")
    ap.add_argument("--matlab", type=float, help="MATLAB time per pair [s]")
    ap.add_argument("--octave", type=float, help="Octave time per pair [s]")
    a = ap.parse_args()
    n = a.size
    mov = ParticleMovie(FLOWS["lorentz_cells"](n, n), (n, n), seed=1)
    fr, _ = mov.run(a.pairs + 1, truth=False)
    res = {}
    if a.v10:
        spec = importlib.util.spec_from_file_location("pex_v10", a.v10)
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        e = m.PExFlow2D()
        e(fr[0], fr[1])
        res["Python v1.0 core"] = timeit(lambda: [e(fr[i], fr[i + 1]) for i in range(a.pairs)], 2) / a.pairs
    cv = pf._cv2
    pf._cv2 = None
    e = pf.PExFlow2D()
    e(fr[0], fr[1])
    res["Python v1.1 (NumPy/SciPy)"] = timeit(lambda: e.sequence(fr), 2) / a.pairs
    pf._cv2 = cv
    if cv is not None:
        e = pf.PExFlow2D()
        e(fr[0], fr[1])
        res["Python v1.1 + OpenCV"] = timeit(lambda: e.sequence(fr), 2) / a.pairs
    if a.matlab:
        res["MATLAB (R2018a+)"] = a.matlab
    if a.octave:
        res["GNU Octave"] = a.octave
    for k, v in res.items():
        print(f"  {k:28s} {v:.3f} s/pair")
    out = os.path.join(ROOT, "docs", "media")
    json.dump(dict(size=n, pairs=a.pairs, cpu=os.cpu_count(), seconds_per_pair=res),
              open(os.path.join(out, "benchmark.json"), "w"), indent=1)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from generate_media import style, BG, FG, MUTED, GRID, ACCENT, ACCENT2
    style()
    names = list(res)
    vals = [res[k] for k in names]
    cols = [MUTED if "v1.0" in k else (ACCENT2 if ("MATLAB" in k or "Octave" in k) else ACCENT)
            for k in names]
    fig, ax = plt.subplots(figsize=(7.6, 0.55 * len(names) + 1.1), dpi=120)
    fig.subplots_adjust(0.3, 0.2, 0.97, 0.8)
    y = np.arange(len(names))[::-1]
    ax.barh(y, vals, color=cols, height=0.62)
    for yi, v in zip(y, vals):
        ax.text(v + max(vals) * 0.01, yi, f"{v:.2f} s", va="center", color=FG, fontsize=9)
    ax.set_yticks(y)
    ax.set_yticklabels(names, color=FG)
    ax.set_xlim(0, max(vals) * 1.15)
    ax.set_xlabel("seconds per image pair (lower is better)")
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(f"Time per pair  ·  {n}x{n} px  ·  defaults  ·  {os.cpu_count()}-thread CPU",
                 fontsize=10.5, loc="right")
    fig.savefig(os.path.join(out, "speed.png"), facecolor=BG)
    print("saved docs/media/speed.png")


if __name__ == "__main__":
    main()
