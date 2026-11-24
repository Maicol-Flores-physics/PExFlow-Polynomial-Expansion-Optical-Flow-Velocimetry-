%% Execution time of PExFlow2D (MATLAB R2018a+ / GNU Octave >= 6)
% Synthetic 512x512 particle sequence, default parameters.
n = 512; T = 11;
if exist('OCTAVE_VERSION', 'builtin'), rand('seed', 1); else, rng(1); end
[X, Y] = meshgrid(0:n - 1, 0:n - 1);
np_ = round(0.06 * n * n);
xp = rand(np_, 1) * n; yp = rand(np_, 1) * n;
k = pi / 64; video = zeros(n, n, T, 'single');
for t = 1:T
    I = zeros(n, n);
    ix = floor(xp); iy = floor(yp);
    for oy = -2:3
        for ox = -2:3
            xx = ix + ox; yy = iy + oy;
            m = xx >= 0 & xx < n & yy >= 0 & yy < n;
            I = I + accumarray([yy(m) + 1, xx(m) + 1], ...
                exp(-((xx(m) - xp(m)) .^ 2 + (yy(m) - yp(m)) .^ 2) / (2 * 0.65 ^ 2)), [n, n]);
        end
    end
    video(:, :, t) = min(230 * I, 255);
    u = 2.5 * sin(k * xp) .* cos(k * yp); v = -2.5 * cos(k * xp) .* sin(k * yp);
    xp = mod(xp + u, n); yp = mod(yp + v, n);
end
PExFlow2D(video(:, :, 1:2));                 % warm-up
tic; PExFlow2D(video); t = toc;
fprintf('PExFlow2D %dx%d: %.3f s per pair (%d pairs)\n', n, n, t / (T - 1), T - 1);
