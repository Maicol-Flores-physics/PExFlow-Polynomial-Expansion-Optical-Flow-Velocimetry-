%% PExFlow 2D -- process a video file and save the velocity fields
% MATLAB R2018a+ (uses VideoReader). For GNU Octave see ../octave.
clear; close all;

%% User parameters -----------------------------------------------------------
video_file  = 'my_experiment.mp4';   % .mp4 / .avi / .mov
dt          = 1 / 30;                % time between consecutive frames [s]
scale       = 1e-4;                  % calibration [m/px]
frame_gap   = 1;                     % K: pairs (t, t+K)
frames_used = [];                    % e.g. 1:200  ([] = all)
roi         = [];                    % [row1 row2 col1 col2] ([] = full frame)

% PExFlow parameters (defaults of the v1.0 application)
window_size = 3;  poly_degree = 2;  pyr_levels = 6;  iterations = 3;  sigma_int = 2.0;

%% Read the video ------------------------------------------------------------
vr = VideoReader(video_file);
frames = {};
while hasFrame(vr)
    f = readFrame(vr);
    if size(f, 3) == 3, f = mean(double(f), 3); end
    frames{end + 1} = single(f); %#ok<SAGROW>
end
video = cat(3, frames{:});
clear frames
if ~isempty(frames_used), video = video(:, :, frames_used); end
if ~isempty(roi),         video = video(roi(1):roi(2), roi(3):roi(4), :); end
fprintf('%d frames of %dx%d px\n', size(video, 3), size(video, 2), size(video, 1));

%% PExFlow -------------------------------------------------------------------
[dX, dY] = PExFlow2D(video, window_size, poly_degree, pyr_levels, iterations, ...
                     sigma_int, frame_gap, 'Verbose', true);
F = pexflow_fields(dX, dY, dt * frame_gap, scale);

%% Save (same fields as the PExFlow application) -----------------------------
U = F.U; V = F.V; x = F.x; y = F.y; t = (0:size(U, 3) - 1) * dt; %#ok<NASGU>
convencion_V = 'positive = upwards (decreasing pixel row)';      %#ok<NASGU>
save('pexflow_result.mat', 'U', 'V', 'x', 'y', 't', 'convencion_V', '-v7.3');

%% Animation: frame | speed + vectors | vorticity ----------------------------
figure('Color', 'w', 'Position', [80 80 1300 420]);
s = 16;
[xs, ys] = meshgrid(1:s:size(video, 2), 1:s:size(video, 1));
wsort = sort(abs(F.vorticity(:)));
wlim = wsort(max(1, round(0.99 * numel(wsort))));    % 99th percentile (no toolbox)
for k = 1:size(U, 3)
    subplot(1, 3, 1); imagesc(video(:, :, k)); axis image off; colormap(gca, gray);
    title(sprintf('frame %d', k));
    subplot(1, 3, 2); imagesc(F.speed(:, :, k)); axis image off; colormap(gca, parula);
    hold on; quiver(xs, ys, U(1:s:end, 1:s:end, k), -V(1:s:end, 1:s:end, k), 'w'); hold off;
    title('|u| [m/s]');
    subplot(1, 3, 3); imagesc(F.vorticity(:, :, k), [-wlim wlim]); axis image off;
    colormap(gca, redblue()); title('\omega_z [1/s]');
    drawnow;
end

function c = redblue()
n = 128; a = linspace(0, 1, n)';
c = [[a, a, ones(n, 1)]; [ones(n, 1), flipud(a), flipud(a)]];
end
