# -*- coding: utf-8 -*-
"""
example_basic.py -- PExFlow 2D in a few lines.

Generates a synthetic particle-image pair of a Lamb-Oseen vortex (known
ground truth), computes the displacement field with PExFlow and plots
velocity, vorticity and the error.

    pip install -r requirements.txt
    python example_basic.py
"""
import os
import sys
import time

import numpy as np
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "examples"))
from pexflow2d import PExFlow2D, to_physical, derived_fields   # noqa: E402
from synthetic import FLOWS, ParticleMovie, r2_score             # noqa: E402

# 1. two particle images with known displacement --------------------------
H, W = 256, 256
movie = ParticleMovie(FLOWS["lamb_oseen"](H, W), shape=(H, W), ppp=0.06, seed=1)
frames, (true_dx, true_dy) = movie.run(n_frames=2)

# 2. PExFlow ----------------------------------------------------------------
est = PExFlow2D(window_size=3, poly_degree=2, pyr_levels=6, iterations=3, sigma_int=2.0)
t0 = time.perf_counter()
dx, dy = est(frames[0], frames[1])                   # pixels, +dy downwards
print(f"PExFlow2D: {time.perf_counter() - t0:.3f} s   "
      f"R^2 = {r2_score(dx, dy, true_dx[0], true_dy[0]):.4f}")

# 3. physical units (dt = 1 ms, 50 um / px) and derived fields ------------
U, V = to_physical(dx, dy, dt=1e-3, scale=50e-6)     # m/s, +V upwards
F = derived_fields(U, V, h=50e-6)

# 4. plot -----------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))
ax[0].imshow(frames[0], cmap="gray")
ax[0].set_title("frame 1")
im = ax[1].imshow(F["speed"], cmap="magma")
s = 12
yy, xx = np.mgrid[s // 2:H:s, s // 2:W:s]
ax[1].quiver(xx, yy, U[yy, xx], V[yy, xx], color="w")
ax[1].set_title("|u| [m/s]")
fig.colorbar(im, ax=ax[1])
im = ax[2].imshow(F["vorticity"], cmap="RdBu_r")
ax[2].set_title(r"$\omega_z$ [1/s]")
fig.colorbar(im, ax=ax[2])
for a in ax:
    a.set_xticks([]); a.set_yticks([])
fig.tight_layout()
plt.show()
