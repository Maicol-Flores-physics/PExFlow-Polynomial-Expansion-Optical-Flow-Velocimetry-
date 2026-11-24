function [dX, dY] = PExFlow2D(video, window_size, poly_degree, pyr_levels, iterations, sigma_int, frame_gap, varargin)
%PEXFLOW2D  Polynomial Expansion Optical Flow Velocimetry (PExFlow 2D).
%
%   [dX, dY] = PExFlow2D(video)
%   [dX, dY] = PExFlow2D(video, window_size, poly_degree, pyr_levels, ...
%                        iterations, sigma_int, frame_gap, 'Name', Value, ...)
%
%   Dense displacement fields between the frames of an image sequence
%   using the PExFlow method (Flores et al. 2026, Exp. Fluids 67:72,
%   https://doi.org/10.1007/s00348-026-04222-1).
%
%   INPUT
%     video        M x N x T array (any numeric class). M x N x 3 x T colour
%                  videos are converted to grey levels.
%     window_size  radius w of the polynomial window, 2w+1 px   (default 3)
%     poly_degree  degree of the local polynomial, >= 2         (default 2)
%                  (odd degrees are promoted to the next even one)
%     pyr_levels   maximum number of pyramid levels             (default 6)
%     iterations   refinements per pyramid level                (default 3)
%     sigma_int    std [px] of the Gaussian integration         (default 2)
%     frame_gap    K: pairs (t, t+K)                            (default 1)
%
%   NAME-VALUE OPTIONS
%     'Reg'        Tikhonov regularisation (relative)            1e-6
%     'Border'     'normalized' (recommended) | 'replicate'      'normalized'
%     'PyrScale'   scale factor between pyramid levels (0,1)     0.65
%     'Presmooth'  OpenCV-style Gaussian pre-smoothing           true
%     'MinSize'    minimum side of the coarsest level [px]       auto
%     'FineBlur'   [0.25 0.5 0.25] blur of the finest level      false
%     'Precision'  'single' (fast) | 'double'                    'single'
%     'Verbose'    print progress                                false
%
%   OUTPUT
%     dX, dY       M x N x (T-K) displacements in PIXELS.
%                  dX: columns (+ to the right). dY: rows (+ DOWNWARDS,
%                  image convention). Physical velocity with +V upwards:
%                      U = dX*scale/dt;   V = -dY*scale/dt;
%
%   Requires MATLAB R2018a or newer, or GNU Octave >= 6 (Octave support is
%   written for compatibility but not yet tested). No toolboxes.
%
%   Optimisations: every frame is expanded once per pyramid level and
%   reused by the two pairs that contain it; border operators, pixel grids
%   and resampling matrices are precomputed once per image size; bilinear
%   warps use precomputed linear indices; separable conv2 everywhere.
%
%   Example
%     [dX, dY] = PExFlow2D(frames, 3, 2, 6, 3, 2.0, 1);
%     quiver(dX(:,:,1), -dY(:,:,1)); axis ij image
%
%   Author: Maicol Flores - Instituto de Fisica, Facultad de Ciencias, UdelaR.
%   See also PEXFLOW_FIELDS.

% ---------------------------------------------------------------- inputs
if nargin < 2 || isempty(window_size), window_size = 3;   end
if nargin < 3 || isempty(poly_degree), poly_degree = 2;   end
if nargin < 4 || isempty(pyr_levels),  pyr_levels  = 6;   end
if nargin < 5 || isempty(iterations),  iterations  = 3;   end
if nargin < 6 || isempty(sigma_int),   sigma_int   = 2.0; end
if nargin < 7 || isempty(frame_gap),   frame_gap   = 1;   end

opt = struct('Reg', 1e-6, 'Border', 'normalized', 'PyrScale', 0.65, ...
             'Presmooth', true, 'MinSize', [], 'FineBlur', false, ...
             'Precision', 'single', 'Verbose', false);
names = fieldnames(opt);
if mod(numel(varargin), 2) ~= 0
    error('PExFlow2D:args', 'Options must be Name-Value pairs.');
end
for k = 1:2:numel(varargin)
    hit = find(strcmpi(names, varargin{k}), 1);
    if isempty(hit)
        error('PExFlow2D:args', 'Unknown option "%s".', varargin{k});
    end
    opt.(names{hit}) = varargin{k + 1};
end

if poly_degree < 2 || poly_degree ~= round(poly_degree)
    error('PExFlow2D:args', 'poly_degree must be an integer >= 2.');
end
if window_size < 1 || window_size ~= round(window_size)
    error('PExFlow2D:args', 'window_size must be an integer >= 1.');
end
if ~(opt.PyrScale > 0 && opt.PyrScale < 1)
    error('PExFlow2D:args', 'PyrScale must be in (0, 1).');
end
border = lower(opt.Border);
if ~any(strcmp(border, {'normalized', 'replicate'}))
    error('PExFlow2D:args', 'Border must be ''normalized'' or ''replicate''.');
end
cls = lower(opt.Precision);
if ~any(strcmp(cls, {'single', 'double'}))
    error('PExFlow2D:args', 'Precision must be ''single'' or ''double''.');
end
if mod(poly_degree, 2) == 1
    poly_degree = poly_degree + 1;                % parity correction
end

if ndims(video) == 4                              % M x N x 3 x T colour
    video = squeeze(mean(double(video), 3));
end
[M, N, T] = size(video);
if T <= frame_gap
    error('PExFlow2D:args', 'The video must have more frames than frame_gap.');
end

% ------------------------------------------------------- precomputation
w  = window_size;
op = build_operator(w, poly_degree, cls);
wk = gauss_kernel(sigma_int, cls);
if isempty(opt.MinSize)
    ms = max(2 * w + 1, 8);
else
    ms = max(opt.MinSize, 2 * w + 1);
end
sizes = pyramid_sizes(M, N, pyr_levels, ms, opt.PyrScale);
nL = size(sizes, 1);
SC = cell(nL, 1);
for l = 1:nL
    SC{l} = shape_cache(op, sizes(l, 1), sizes(l, 2), border, cls);
    if l < nL      % resampling matrices: level l -> level l+1 (flow)
        SC{l}.Uy = interp_matrix(sizes(l, 1), sizes(l + 1, 1));
        SC{l}.Ux = interp_matrix(sizes(l, 2), sizes(l + 1, 2));
    end
    if opt.Presmooth   % original -> level l (pyramid of the images)
        SC{l}.Dy = interp_matrix(M, sizes(l, 1));
        SC{l}.Dx = interp_matrix(N, sizes(l, 2));
    else
        SC{l}.Dy = area_matrix(M, sizes(l, 1));
        SC{l}.Dx = area_matrix(N, sizes(l, 2));
    end
end
P = struct('op', op, 'wk', wk, 'iter', max(1, iterations), 'reg', opt.Reg, ...
           'border', border, 'presmooth', logical(opt.Presmooth), ...
           'fineblur', logical(opt.FineBlur), 'cls', cls, 'sizes', sizes, 'N', N);

% ------------------------------------------------------------ main loop
n  = T - frame_gap;
dX = zeros(M, N, n, cls);
dY = zeros(M, N, n, cls);
E  = cell(T, 1);                                  % expansion cache
t0 = tic;
for t = 1:n
    for kk = [t, t + frame_gap]
        if isempty(E{kk})
            E{kk} = expand_frame(cast(video(:, :, kk), cls), P, SC);
        end
    end
    [dX(:, :, t), dY(:, :, t)] = solve_pair(E{t}, E{t + frame_gap}, P, SC);
    E{t} = [];
    if opt.Verbose
        el = toc(t0);
        fprintf('  PExFlow2D  pair %d/%d  (%.1f s, %.2f s/pair)\n', t, n, el, el / t);
    end
end
end

% =======================================================================
%  Operator
% =======================================================================
function op = build_operator(w, deg, cls)
ya = []; xb = [];
for d = 0:deg
    for a = d:-1:0
        ya(end + 1) = a;      %#ok<AGROW>
        xb(end + 1) = d - a;  %#ok<AGROW>
    end
end
nT  = numel(ya);
r   = (-w:w)';
u   = r / w;
app = exp(-r .^ 2 / (2 * (0.3 * w) ^ 2));

eqs = zeros(0, 5); used = [];
for i = 1:nT
    a = ya(i); b = xb(i);
    if a + b < 1 || a + b > deg - 1, continue; end
    iy = find(ya == a + 1 & xb == b);
    ix = find(ya == a & xb == b + 1);
    eqs(end + 1, :) = [i, iy, a + 1, ix, b + 1];  %#ok<AGROW>
    used = [used, i, iy, ix];                     %#ok<AGROW>
end
k = zeros(2 * w + 1, deg + 1);
for p = 0:deg
    k(:, p + 1) = app .* u .^ p;
end
op = struct('w', w, 'deg', deg, 'ya', ya, 'xb', xb, 'nT', nT, 'u', u, ...
            'app', app, 'eqs', eqs, 'used', unique(used), 'cls', cls, ...
            'k', cast(k, cls), 'scale', cast(1 ./ w .^ (ya + xb), cls));
op.C = class_operator(op, -w, w, -w, w);
end

function C = class_operator(op, rlo, rhi, clo, chi)
w  = op.w;
t  = (-w:w)';
p2 = 0:(2 * op.deg);
wy = op.app .* (t >= rlo & t <= rhi);
wx = op.app .* (t >= clo & t <= chi);
my = sum(bsxfun(@times, wy, bsxfun(@power, op.u, p2)), 1);
mx = sum(bsxfun(@times, wx, bsxfun(@power, op.u, p2)), 1);
G  = my(bsxfun(@plus, op.ya', op.ya) + 1) .* mx(bsxfun(@plus, op.xb', op.xb) + 1);
C  = cast(pinv(G), op.cls);
end

function sc = shape_cache(op, H, W, border, cls)
w = op.w;
[rlo, rhi] = border_classes(H, w);
[clo, chi] = border_classes(W, w);
sc.H = H; sc.W = W;
[X, Y] = meshgrid(0:W - 1, 0:H - 1);
sc.X = cast(X, cls); sc.Y = cast(Y, cls);
if ~strcmp(border, 'normalized'), return; end
sc.rb = find(rlo > -w | rhi < w);
sc.cb = find(clo > -w | chi < w);
sc.rows_in = setdiff(1:H, sc.rb);
[ucls, ~, cidx] = unique([clo(:), chi(:)], 'rows');
ncls = size(ucls, 1);
sc.colsets = cell(ncls, 1);
for j = 1:ncls
    sc.colsets{j} = find(cidx == j);
end
sc.row_ops = cell(numel(sc.rb), ncls);
for ii = 1:numel(sc.rb)
    r = sc.rb(ii);
    for j = 1:ncls
        sc.row_ops{ii, j} = class_operator(op, rlo(r), rhi(r), ucls(j, 1), ucls(j, 2)).';
    end
end
sc.col_ops = cell(numel(sc.cb), 1);
for jj = 1:numel(sc.cb)
    c = sc.cb(jj);
    sc.col_ops{jj} = class_operator(op, -w, w, clo(c), chi(c)).';
end
end

function [lo, hi] = border_classes(L, w)
idx = 0:L - 1;
lo = max(-w, -idx);
hi = min(w, L - 1 - idx);
end

% =======================================================================
%  Pyramid helpers
% =======================================================================
function sizes = pyramid_sizes(M, N, L, min_size, scale)
sizes = [M, N];
s = 1;
for k = 2:L
    s = s * scale;
    sz = [round_even(M * s), round_even(N * s)];
    if min(sz) < min_size || isequal(sz, sizes(1, :)), break; end
    sizes = [sz; sizes]; %#ok<AGROW>
end
end

function r = round_even(x)          % round half to even (as Python)
r = round(x);
if abs(x - fix(x)) == 0.5
    r = 2 * round(x / 2);
end
end

function R = interp_matrix(n_in, n_out)
% Bilinear resampling matrix (pixel-centre convention, clamped borders).
c = ((0:n_out - 1)' + 0.5) * n_in / n_out - 0.5;
c = min(max(c, 0), n_in - 1);
if n_in == 1
    R = sparse(ones(n_out, 1)); return;
end
i0 = min(floor(c), n_in - 2);
f  = c - i0;
R  = sparse([(1:n_out)'; (1:n_out)'], [i0 + 1; i0 + 2], [1 - f; f], n_out, n_in);
end

function R = area_matrix(n_in, n_out)
sc = n_in / n_out;
I = []; J = []; V = [];
for d = 0:n_out - 1
    a = d * sc; b = (d + 1) * sc;
    ii = floor(a):min(ceil(b), n_in) - 1;
    v  = min(b, ii + 1) - max(a, ii);
    I = [I; repmat(d + 1, numel(ii), 1)]; %#ok<AGROW>
    J = [J; ii(:) + 1];                   %#ok<AGROW>
    V = [V; v(:) / sum(v)];               %#ok<AGROW>
end
R = sparse(I, J, V, n_out, n_in);
end

function out = resample(img, Ry, Rx, cls)
out = cast(Ry * double(img) * Rx.', cls);
end

function k = gauss_kernel(sigma, cls)
r = max(1, ceil(3 * sigma));
x = (-r:r)';
k = exp(-x .^ 2 / (2 * sigma ^ 2));
k = cast(k / sum(k), cls);
end

% =======================================================================
%  Filtering (correlation, SciPy/OpenCV semantics)
% =======================================================================
function idx = pad_index(n, r, mode)
i = (1 - r):(n + r);
if strcmp(mode, 'replicate')
    idx = min(max(i, 1), n);
else                                              % 'mirror' (reflect-101)
    p = 2 * (n - 1);
    j = mod(i - 1, p);
    j(j >= n) = p - j(j >= n);
    idx = j + 1;
end
end

function out = corr1(img, k, dim, mode)
% 1-D correlation of img along dim with kernel k ('zero'|'replicate'|'mirror').
kf = flipud(k(:));
r  = (numel(kf) - 1) / 2;
if dim == 1
    if strcmp(mode, 'zero')
        out = conv2(kf, 1, img, 'same');
    else
        out = conv2(kf, 1, img(pad_index(size(img, 1), r, mode), :), 'valid');
    end
else
    if strcmp(mode, 'zero')
        out = conv2(1, kf.', img, 'same');
    else
        out = conv2(1, kf.', img(:, pad_index(size(img, 2), r, mode)), 'valid');
    end
end
end

function out = smooth2(img, k, mode)
kf = flipud(k(:));
r  = (numel(kf) - 1) / 2;
if strcmp(mode, 'zero')
    out = conv2(kf, kf.', img, 'same');
else
    out = conv2(kf, kf.', img(pad_index(size(img, 1), r, mode), ...
                              pad_index(size(img, 2), r, mode)), 'valid');
end
end

function img = gaussian_presmooth(img, sigma, cls)
if sigma <= 0, return; end
ks = round(sigma * 5);
if mod(ks, 2) == 0, ks = ks + 1; end
ks = max(ks, 3);
r  = (ks - 1) / 2;
x  = (-r:r)';
k  = exp(-x .^ 2 / (2 * sigma ^ 2));
k  = cast(k / sum(k), cls);
img = smooth2(img, k, 'mirror');
end

% =======================================================================
%  Polynomial expansion of one frame (all pyramid levels)
% =======================================================================
function E = expand_frame(f, P, SC)
if P.fineblur
    f = smooth2(f, cast([0.25; 0.5; 0.25], P.cls), 'mirror');
end
nL = size(P.sizes, 1);
E  = cell(nL, 1);
for l = 1:nL
    if l == nL
        I = f;
    elseif P.presmooth
        s = P.sizes(l, 2) / P.N;
        I = resample(gaussian_presmooth(f, (1 / s - 1) * 0.5, P.cls), SC{l}.Dy, SC{l}.Dx, P.cls);
    else
        I = resample(f, SC{l}.Dy, SC{l}.Dx, P.cls);
    end
    E{l} = poly_expand(I, P.op, SC{l}, P.border);
end
end

function coef = poly_expand(img, op, sc, border)
[H, W] = size(img);
nT = op.nT;
if strcmp(border, 'normalized'), mode = 'zero'; else, mode = 'replicate'; end
Vy = cell(op.deg + 1, 1);
for a = 0:op.deg
    Vy{a + 1} = corr1(img, op.k(:, a + 1), 1, mode);
end
Pm = zeros(H, W, nT, op.cls);
for i = 1:nT
    Pm(:, :, i) = corr1(Vy{op.ya(i) + 1}, op.k(:, op.xb(i) + 1), 2, mode);
end
coef = reshape(reshape(Pm, H * W, nT) * op.C.', H, W, nT);
if strcmp(border, 'normalized')
    for ii = 1:numel(sc.rb)                       % border rows
        r  = sc.rb(ii);
        Pr = reshape(Pm(r, :, :), W, nT);
        for j = 1:numel(sc.colsets)
            cols = sc.colsets{j};
            coef(r, cols, :) = reshape(Pr(cols, :) * sc.row_ops{ii, j}, 1, numel(cols), nT);
        end
    end
    ri = sc.rows_in;
    for jj = 1:numel(sc.cb)                       % border columns
        c  = sc.cb(jj);
        Pc = reshape(Pm(ri, c, :), numel(ri), nT);
        coef(ri, c, :) = reshape(Pc * sc.col_ops{jj}, numel(ri), 1, nT);
    end
end
coef = bsxfun(@times, coef, reshape(op.scale, 1, 1, nT));
end

% =======================================================================
%  Coarse-to-fine solver
% =======================================================================
function [dx, dy] = solve_pair(E1, E2, P, SC)
nL = size(P.sizes, 1);
dx = zeros(P.sizes(1, :), P.cls);
dy = zeros(P.sizes(1, :), P.cls);
for l = 1:nL
    if l > 1
        sy = P.sizes(l, 1) / P.sizes(l - 1, 1);
        sx = P.sizes(l, 2) / P.sizes(l - 1, 2);
        dx = resample(dx, SC{l - 1}.Uy, SC{l - 1}.Ux, P.cls) * sx;
        dy = resample(dy, SC{l - 1}.Uy, SC{l - 1}.Ux, P.cls) * sy;
    end
    for it = 1:P.iter
        [dx, dy] = solve_level(E1{l}, E2{l}, dx, dy, P, SC{l});
    end
end
end

function [dx, dy] = solve_level(c1, c2, dxp, dyp, P, sc)
op = P.op; cls = P.cls;
[H, W] = size(dxp);
yb = sc.Y + dyp;
xb = sc.X + dxp;

% bilinear warp of the coefficient planes (indices computed once)
yc = min(max(yb, 0), H - 1);
xc = min(max(xb, 0), W - 1);
y0 = min(floor(yc), max(H - 2, 0));
x0 = min(floor(xc), max(W - 2, 0));
wy = yc - y0;
wx = xc - x0;
i00 = double(y0) + double(x0) * H + 1;
i01 = i00 + H * (W > 1);
i10 = i00 + (H > 1);
i11 = i10 + H * (W > 1);
w11 = wy .* wx;  w10 = wy - w11;  w01 = wx - w11;  w00 = 1 - wy - w01;
c2w = zeros(H, W, op.nT, cls);
for p = op.used
    q = c2(:, :, p);
    c2w(:, :, p) = w00 .* q(i00) + w01 .* q(i01) + w10 .* q(i10) + w11 .* q(i11);
end

normalized = strcmp(P.border, 'normalized');
if normalized
    valid = cast(xb >= 0 & xb <= W - 1 & yb >= 0 & yb <= H - 1, cls);
    R   = op.w + 1;
    dxe = min(min(sc.X, xb), min(W - 1 - sc.X, W - 1 - xb));
    dye = min(min(sc.Y, yb), min(H - 1 - sc.Y, H - 1 - yb));
    rx  = min(max((dxe + 1) / R, 0), 1);
    ry  = min(max((dye + 1) / R, 0), 1);
    cw  = valid .* (rx .* ry) .^ 4;
    smode = 'zero';
else
    smode = 'replicate';
end

G11 = zeros(H, W, cls); G12 = G11; G22 = G11; h1 = G11; h2 = G11;
for e = 1:size(op.eqs, 1)
    i0 = op.eqs(e, 1); iy = op.eqs(e, 2); fy = op.eqs(e, 3);
    ix = op.eqs(e, 4); fx = op.eqs(e, 5);
    Ay = (0.5 * fy) * (c1(:, :, iy) + c2w(:, :, iy));
    Ax = (0.5 * fx) * (c1(:, :, ix) + c2w(:, :, ix));
    Db = c1(:, :, i0) - c2w(:, :, i0) + Ay .* dyp + Ax .* dxp;
    G11 = G11 + Ay .* Ay;  G12 = G12 + Ay .* Ax;  G22 = G22 + Ax .* Ax;
    h1  = h1  + Ay .* Db;  h2  = h2  + Ax .* Db;
end
if normalized
    G11 = G11 .* cw; G12 = G12 .* cw; G22 = G22 .* cw; h1 = h1 .* cw; h2 = h2 .* cw;
end
G11 = smooth2(G11, P.wk, smode); G12 = smooth2(G12, P.wk, smode);
G22 = smooth2(G22, P.wk, smode);
h1  = smooth2(h1,  P.wk, smode); h2  = smooth2(h2,  P.wk, smode);

lam = P.reg * mean(G11(:) + G22(:)) + realmin(cls);
G11 = G11 + lam; G22 = G22 + lam;
dt  = G11 .* G22 - G12 .^ 2;
dy  = (G22 .* h1 - G12 .* h2) ./ dt;
dx  = (G11 .* h2 - G12 .* h1) ./ dt;

if normalized
    alpha = min(max(smooth2(valid, P.wk, 'replicate'), 0), 1);
    if min(alpha(:)) < 1
        kb  = gauss_kernel(4.0, cls);
        den = smooth2(alpha, kb, 'zero') + 1e-12;
        fxv = smooth2(alpha .* dx, kb, 'zero') ./ den;
        fyv = smooth2(alpha .* dy, kb, 'zero') ./ den;
        dx  = alpha .* dx + (1 - alpha) .* fxv;
        dy  = alpha .* dy + (1 - alpha) .* fyv;
    end
end
end
