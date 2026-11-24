# -*- coding: utf-8 -*-
"""Tests for pexflow2d.py   --   run with:  python -m pytest python/tests
(or simply  python python/tests/test_pexflow2d.py)."""
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "python"))
sys.path.insert(0, os.path.join(ROOT, "examples"))

from pexflow2d import PExFlow2D, pexflow2d, to_physical, derived_fields, _main  # noqa: E402
from synthetic import FLOWS, ParticleMovie, r2_score                           # noqa: E402


def _pair(flow, shape=(160, 192), seed=0, **kw):
    mov = ParticleMovie(FLOWS[flow](*shape, **kw), shape, seed=seed)
    fr, (tx, ty) = mov.run(3)
    return fr, tx, ty


def test_uniform_translation():
    def vel(x, y, t):
        return 1.7 + 0 * x, -0.6 + 0 * y
    mov = ParticleMovie(vel, (128, 128), seed=1)
    fr, _ = mov.run(2, truth=False)
    dx, dy = PExFlow2D()(fr[0], fr[1])
    c = (slice(8, -8), slice(8, -8))
    assert abs(np.median(dx[c]) - 1.7) < 0.03
    assert abs(np.median(dy[c]) + 0.6) < 0.03


def test_lamb_oseen_accuracy():
    fr, tx, ty = _pair("lamb_oseen")
    dx, dy = PExFlow2D()(fr[0], fr[1])
    assert r2_score(dx, dy, tx[0], ty[0], 4) > 0.99


def test_parameters_and_degree_parity():
    fr, tx, ty = _pair("lorentz_cells", seed=2)
    for w, d in [(2, 3), (4, 4), (5, 6)]:
        est = PExFlow2D(window_size=w, poly_degree=d)
        assert est.params["poly_degree"] % 2 == 0
        dx, dy = est(fr[0], fr[1])
        assert r2_score(dx, dy, tx[0], ty[0], 4) > 0.97


def test_sequence_equals_pairs_and_matlab_layout():
    fr, _, _ = _pair("double_gyre", seed=3)
    est = PExFlow2D()
    dX, dY = est.sequence(fr)
    dx1, dy1 = est(fr[1], fr[2])
    assert np.array_equal(dX[1], dx1) and np.array_equal(dY[1], dy1)
    mX, mY = pexflow2d(np.moveaxis(fr, 0, 2))
    assert mX.shape == (fr.shape[1], fr.shape[2], 2)
    assert np.allclose(mX[:, :, 1], dx1)


def test_replicate_border_runs():
    fr, tx, ty = _pair("poiseuille", seed=4)
    dx, dy = PExFlow2D(border="replicate", precision="float64")(fr[0], fr[1])
    assert dx.dtype == np.float64
    assert r2_score(dx, dy, tx[0], ty[0], 8) > 0.98


def test_derived_fields_rigid_rotation():
    H = W = 64
    Y, X = np.mgrid[0:H, 0:W].astype(float)
    om = 0.02                                   # physical: x right, y up = -row
    U = -om * (-(Y - H / 2))
    V = om * (X - W / 2)
    F = derived_fields(U, V)
    assert np.allclose(F["vorticity"], 2 * om)
    assert np.allclose(F["divergence"], 0, atol=1e-12)
    u, v = to_physical(np.ones((2, 2)), np.ones((2, 2)), dt=0.5, scale=2.0)
    assert np.all(u == 4) and np.all(v == -4)


def test_invalid_arguments():
    for kw in [dict(poly_degree=1), dict(window_size=0), dict(pyr_scale=1.2),
               dict(border="wrap"), dict(precision="float16")]:
        try:
            PExFlow2D(**kw)
        except ValueError:
            continue
        raise AssertionError(f"{kw} should raise ValueError")


def test_cli_roundtrip():
    import imageio.v3 as iio
    fr, _, _ = _pair("lamb_oseen", shape=(96, 96), seed=5)
    with tempfile.TemporaryDirectory() as d:
        for i, f in enumerate(fr):
            iio.imwrite(os.path.join(d, f"f{i:03d}.png"), f.astype(np.uint8))
        out = os.path.join(d, "res.npz")
        _main([os.path.join(d, "*.png"), "--out", out, "--dt", "0.5", "--scale", "1e-3"])
        with np.load(out) as r:          # close the file (Windows)
            assert r["U"].shape == (96, 96, 2) and r["x"].size == 96


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok ", name)
