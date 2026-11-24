# -*- coding: utf-8 -*-
"""
synthetic.py -- synthetic particle-image generator with known ground truth.

Used by the examples, the validation tests and to produce the media of the
README. Each flow is a function  vel(x, y, t) -> (u, v)  in PIXEL units per
frame, with image coordinates (x to the right, y DOWNWARDS).

    from synthetic import FLOWS, ParticleMovie
    mov = ParticleMovie(FLOWS["lamb_oseen"], shape=(256, 256), ppp=0.06)
    frames, (dx_true, dy_true) = mov.run(n_frames=20)
"""
import numpy as np


# ---------------------------------------------------------------------------
# Canonical flows  (pixel / frame, image coordinates)
# ---------------------------------------------------------------------------
def lamb_oseen(H, W, gamma=2200.0, rc=28.0):
    xc, yc = W / 2, H / 2

    def vel(x, y, t):
        dx, dy = x - xc, y - yc
        r2 = dx ** 2 + dy ** 2 + 1e-9
        ut = gamma / (2 * np.pi * r2) * (1 - np.exp(-r2 / rc ** 2))
        return -ut * dy, ut * dx
    return vel


def lorentz_cells(H, W, cell=64.0, U0=2.5, period=60.0):
    """Chessboard array of counter-rotating vortices (Taylor-Green-like),
    similar to the Lorentz-forced experiment of the paper, with a slow
    oscillation of the amplitude."""
    k = np.pi / cell

    def vel(x, y, t):
        a = U0 * (0.75 + 0.25 * np.cos(2 * np.pi * t / period))
        u = a * np.sin(k * x) * np.cos(k * y)
        v = -a * np.cos(k * x) * np.sin(k * y)
        return u, v
    return vel


def poiseuille(H, W, umax=4.0):
    def vel(x, y, t):
        s = y / (H - 1)
        return umax * 4 * s * (1 - s), 0 * x
    return vel


def sink(H, W, q=700.0, gamma=700.0, core=200.0):
    """Sink + swirl with a regularised core (max |u| ~ 5.5 px/frame)."""
    xc, yc = W / 2, H / 2

    def vel(x, y, t):
        dx, dy = x - xc, y - yc
        r2 = dx ** 2 + dy ** 2 + core
        return (-q * dx - gamma * dy) / (2 * np.pi * r2), (-q * dy + gamma * dx) / (2 * np.pi * r2)
    return vel


def turbulence(H, W, slope=-3.0, kmin=3, rms=1.6, seed=1):
    """Frozen divergence-free random field with E(k) ~ k^slope
    (2-D 'Kraichnan-like' synthetic turbulence) -- streamfunction based."""
    rng = np.random.default_rng(seed)
    ky = np.fft.fftfreq(H) * H              # isotropic: cycles per image height
    kx = np.fft.fftfreq(W) * H
    KX, KY = np.meshgrid(kx, ky)
    K = np.hypot(KX, KY)
    K[0, 0] = 1
    amp = np.where(K >= kmin, K ** ((slope - 1) / 2) / K, 0.0)   # |psi_k|
    psi_k = amp * np.exp(2j * np.pi * rng.random((H, W)))
    psi = np.real(np.fft.ifft2(psi_k))
    uu = np.gradient(psi, axis=0)          # u =  dpsi/dy
    vv = -np.gradient(psi, axis=1)         # v = -dpsi/dx
    s = rms / np.sqrt(np.mean(uu ** 2 + vv ** 2))
    uu *= s
    vv *= s
    from scipy.ndimage import map_coordinates

    def vel(x, y, t):
        c = [np.mod(y, H), np.mod(x, W)]
        return (map_coordinates(uu, c, order=3, mode="grid-wrap"),
                map_coordinates(vv, c, order=3, mode="grid-wrap"))
    vel.field = (uu, vv)
    return vel


def double_gyre(H, W, eps=0.25, period=40.0, speed=3.0):
    """Time-periodic double gyre on [0,2]x[0,1] mapped to the image
    (max speed ~ ``speed`` px/frame)."""
    om = 2 * np.pi / period

    def vel(x, y, t):
        X = 2.0 * x / W
        Y = 1.0 - y / H
        a = eps * np.sin(om * t)
        b = 1 - 2 * eps * np.sin(om * t)
        f = a * X ** 2 + b * X
        dfx = 2 * a * X + b
        u = -np.pi * np.sin(np.pi * f) * np.cos(np.pi * Y)
        v = np.pi * np.cos(np.pi * f) * np.sin(np.pi * Y) * dfx
        return speed / np.pi * u, -speed / np.pi * v       # y image = downwards
    return vel


def vortex_street(H, W, U=2.2, a=90.0, gamma=520.0, rc=11.0, uc=0.8):
    """Kinematic von Karman-type vortex street: uniform stream U plus two
    staggered rows of Lamb-Oseen vortices of opposite sign convected at
    uc*U. Vortices fade in/out smoothly outside the field of view."""
    h = 0.28 * a * 1.4
    m = 4 * a                                             # chain extends 4a beyond the view
    L = int(np.ceil((W + 2 * m) / a)) * a                 # periodic span (multiple of a)
    ks = np.arange(int(L / a))

    def vel(x, y, t):
        u = np.full(np.shape(x), U, float)
        v = np.zeros(np.shape(x), float)
        for row, sgn in ((0, 1.0), (1, -1.0)):
            yc = H / 2 + (h / 2 if row == 0 else -h / 2)
            for k in ks:
                xc = np.mod(k * a + row * a / 2 + uc * U * t, L) - m
                fade = np.clip((xc + m) / a, 0, 1) * np.clip((L - m - xc) / a, 0, 1)
                fade = fade * fade * (3 - 2 * fade)
                dx, dy = x - xc, y - yc
                r2 = dx ** 2 + dy ** 2 + 1e-9
                ut = sgn * fade * gamma / (2 * np.pi * r2) * (1 - np.exp(-r2 / rc ** 2))
                u = u - ut * dy
                v = v + ut * dx
        return u, v
    return vel


FLOWS = {
    "vortex_street": vortex_street,
    "lamb_oseen": lamb_oseen,
    "double_gyre": double_gyre,
    "lorentz_cells": lorentz_cells,
    "poiseuille": poiseuille,
    "sink": sink,
    "turbulence": turbulence,
}


# ---------------------------------------------------------------------------
# Particle renderer
# ---------------------------------------------------------------------------
def render(px, py, H, W, diam=2.6, intensity=None, noise=0.0, rng=None):
    """Gaussian particle images (float32 in [0, 255])."""
    img = np.zeros((H, W), np.float64)
    s = diam / 4.0
    r = int(np.ceil(3 * s)) + 1
    ix = np.floor(px).astype(int)
    iy = np.floor(py).astype(int)
    I0 = np.ones_like(px) if intensity is None else intensity
    for oy in range(-r, r + 2):
        for ox in range(-r, r + 2):
            xx, yy = ix + ox, iy + oy
            m = (xx >= 0) & (xx < W) & (yy >= 0) & (yy < H)
            val = I0[m] * np.exp(-((xx[m] - px[m]) ** 2 + (yy[m] - py[m]) ** 2) / (2 * s * s))
            np.add.at(img, (yy[m], xx[m]), val)
    img = 230.0 * img
    if noise > 0:
        rng = rng or np.random.default_rng()
        img = img + rng.normal(0, noise, img.shape)
    return np.clip(img, 0, 255).astype(np.float32)


def rk4(vel, x, y, t, h=1.0, n=4):
    dt = h / n
    for i in range(n):
        tt = t + i * dt
        k1 = vel(x, y, tt)
        k2 = vel(x + 0.5 * dt * k1[0], y + 0.5 * dt * k1[1], tt + 0.5 * dt)
        k3 = vel(x + 0.5 * dt * k2[0], y + 0.5 * dt * k2[1], tt + 0.5 * dt)
        k4 = vel(x + dt * k3[0], y + dt * k3[1], tt + dt)
        x = x + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        y = y + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
    return x, y


class ParticleMovie:
    """Lagrangian particle movie. Particles leaving the domain are re-seeded
    at random positions (as in a real light sheet)."""

    def __init__(self, vel, shape=(256, 256), ppp=0.06, diam=2.6, noise=2.0,
                 periodic=False, inflow=None, seed=0):
        self.vel, self.H, self.W = vel, *shape
        self.diam, self.noise, self.periodic = diam, noise, periodic
        self.inflow = inflow                # None | 'left': re-seed in a strip at x < 0
        self.rng = np.random.default_rng(seed)
        n = int(ppp * self.H * self.W)
        m = 6
        self.x = self.rng.uniform(-m, self.W + m, n)
        self.y = self.rng.uniform(-m, self.H + m, n)
        self.I = self.rng.uniform(0.6, 1.0, n)

    def frame(self):
        return render(self.x, self.y, self.H, self.W, self.diam, self.I, self.noise, self.rng)

    def true_displacement(self, t):
        """Exact displacement over one frame of the pixel grid (RK4)."""
        Y, X = np.mgrid[0:self.H, 0:self.W].astype(float)
        x1, y1 = rk4(self.vel, X, Y, t, 1.0, 8)
        return (x1 - X).astype(np.float32), (y1 - Y).astype(np.float32)

    def step(self, t):
        self.x, self.y = rk4(self.vel, self.x, self.y, t, 1.0, 4)
        m = 6
        if self.periodic:
            self.x = np.mod(self.x + m, self.W + 2 * m) - m
            self.y = np.mod(self.y + m, self.H + 2 * m) - m
        else:
            out = (self.x < -m) | (self.x > self.W + m) | (self.y < -m) | (self.y > self.H + m)
            k = int(out.sum())
            if k and self.inflow == "left":
                # re-enter upstream, keeping the Lagrangian spacing (x - span)
                self.x[out & (self.x > self.W + m)] -= self.W + 2 * m
                rest = out & ((self.x < -m - 1e-9) | (self.y < -m) | (self.y > self.H + m))
                j = int(rest.sum())
                self.x[rest] = self.rng.uniform(-m, 0, j)
                self.y[rest] = self.rng.uniform(-m, self.H + m, j)
            elif k:
                self.x[out] = self.rng.uniform(-m, self.W + m, k)
                self.y[out] = self.rng.uniform(-m, self.H + m, k)

    def run(self, n_frames=2, truth=True):
        frames, tdx, tdy = [], [], []
        for t in range(n_frames):
            frames.append(self.frame())
            if truth and t < n_frames - 1:
                a, b = self.true_displacement(t)
                tdx.append(a)
                tdy.append(b)
            self.step(t)
        return np.array(frames), (np.array(tdx), np.array(tdy))


def r2_score(dx, dy, tx, ty, crop=0):
    """2-D coefficient of determination used in the paper."""
    if crop:
        sl = (slice(crop, -crop), slice(crop, -crop))
        dx, dy, tx, ty = dx[sl], dy[sl], tx[sl], ty[sl]
    ss_res = np.sum((dx - tx) ** 2 + (dy - ty) ** 2)
    ss_tot = np.sum((tx - tx.mean()) ** 2 + (ty - ty.mean()) ** 2)
    return 1 - ss_res / ss_tot


def cosine_similarity(dx, dy, tx, ty, crop=0):
    if crop:
        sl = (slice(crop, -crop), slice(crop, -crop))
        dx, dy, tx, ty = dx[sl], dy[sl], tx[sl], ty[sl]
    num = dx * tx + dy * ty
    den = np.hypot(dx, dy) * np.hypot(tx, ty) + 1e-12
    return float(np.mean(num / den))
