%% PExFlow 2D -- synthetic example for GNU Octave (>= 6)
% Run from this folder:   octave --no-gui example_synthetic.m
% Generates a particle-image pair of a Lamb-Oseen vortex with KNOWN
% displacement, runs PExFlow2D and compares with the ground truth.
1;   % (Octave: this file is a script, the functions below are local)

function [u, v] = lamb_oseen(x, y, xc, yc, G, rc)
  dx = x - xc; dy = y - yc;
  r2 = dx .^ 2 + dy .^ 2 + 1e-9;
  ut = G ./ (2 * pi * r2) .* (1 - exp(-r2 / rc ^ 2));
  u = -ut .* dy;  v = ut .* dx;
end

function [x, y] = rk4(vel, x, y, h)
  n = 8; d = h / n;
  for i = 1:n
    [k1u, k1v] = vel(x, y);
    [k2u, k2v] = vel(x + 0.5 * d * k1u, y + 0.5 * d * k1v);
    [k3u, k3v] = vel(x + 0.5 * d * k2u, y + 0.5 * d * k2v);
    [k4u, k4v] = vel(x + d * k3u, y + d * k3v);
    x = x + d / 6 * (k1u + 2 * k2u + 2 * k3u + k4u);
    y = y + d / 6 * (k1v + 2 * k2v + 2 * k3v + k4v);
  end
end

function I = render_particles(x, y, I0, H, W, diam)
  s = diam / 4; r = ceil(3 * s) + 1;
  I = zeros(H, W);
  ix = floor(x); iy = floor(y);
  for oy = -r:r + 1
    for ox = -r:r + 1
      xx = ix + ox; yy = iy + oy;
      m = xx >= 0 & xx < W & yy >= 0 & yy < H;
      val = I0(m) .* exp(-((xx(m) - x(m)) .^ 2 + (yy(m) - y(m)) .^ 2) / (2 * s ^ 2));
      I = I + accumarray([yy(m) + 1, xx(m) + 1], val, [H, W]);
    end
  end
  I = max(min(230 * I + 2 * randn(H, W), 255), 0);
end

%% 1. Synthetic particle images ------------------------------------------
H = 256; W = 256; ppp = 0.06; diam = 2.6; Gamma = 2200; rc = 28;
rand('seed', 1); randn('seed', 1);
np_ = round(ppp * H * W);
xp = rand(np_, 1) * (W + 12) - 6;
yp = rand(np_, 1) * (H + 12) - 6;
Ip = 0.6 + 0.4 * rand(np_, 1);
vel = @(x, y) lamb_oseen(x, y, W / 2, H / 2, Gamma, rc);
[xq, yq] = rk4(vel, xp, yp, 1);
I1 = render_particles(xp, yp, Ip, H, W, diam);
I2 = render_particles(xq, yq, Ip, H, W, diam);
[X, Y] = meshgrid(0:W - 1, 0:H - 1);
[X1, Y1] = rk4(vel, X, Y, 1);
tx = X1 - X;  ty = Y1 - Y;

%% 2. PExFlow --------------------------------------------------------------
tic;
[dX, dY] = PExFlow2D(cat(3, I1, I2), 3, 2, 6, 3, 2.0, 1);
printf('PExFlow2D: %.3f s for one %dx%d pair\n', toc, H, W);

%% 3. Accuracy --------------------------------------------------------------
ex = double(dX) - tx;  ey = double(dY) - ty;
R2 = 1 - sum(ex(:) .^ 2 + ey(:) .^ 2) / ...
         sum((tx(:) - mean(tx(:))) .^ 2 + (ty(:) - mean(ty(:))) .^ 2);
printf('R^2 = %.4f   RMS error = %.3f px\n', R2, sqrt(mean(ex(:) .^ 2 + ey(:) .^ 2)));

%% 4. Figure ----------------------------------------------------------------
F = pexflow_fields(dX, dY, 1, 1);
figure('Color', 'w');
subplot(1, 2, 1); imagesc(I1); axis image; colormap(gray); title('Frame 1');
subplot(1, 2, 2); imagesc(F.speed); axis image; colorbar; hold on;
s = 12; [xs, ys] = meshgrid(1:s:W, 1:s:H);
quiver(xs, ys, dX(1:s:end, 1:s:end), dY(1:s:end, 1:s:end), 1.5, 'w');
title(sprintf('PExFlow |u| [px/frame]   R^2 = %.4f', R2));
