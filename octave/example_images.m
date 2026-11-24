%% PExFlow 2D -- process a folder of images with GNU Octave (>= 6)
% Export your video as an image sequence first, e.g. with ffmpeg:
%     ffmpeg -i my_experiment.mp4 frames/frame_%05d.png
% Then run:   octave --no-gui example_images.m
1;

%% User parameters -----------------------------------------------------------
pattern   = fullfile('frames', '*.png');   % image sequence (sorted by name)
dt        = 1 / 30;                        % time between frames [s]
scale     = 1e-4;                          % calibration [m/px]
frame_gap = 1;                             % K: pairs (t, t+K)
window_size = 3; poly_degree = 2; pyr_levels = 6; iterations = 3; sigma_int = 2.0;

%% Read ------------------------------------------------------------------------
files = dir(pattern);
if isempty(files), error('No images match %s', pattern); end
[~, order] = sort({files.name});
files = files(order);
I = imread(fullfile(files(1).folder, files(1).name));
video = zeros(size(I, 1), size(I, 2), numel(files), 'single');
for k = 1:numel(files)
  I = double(imread(fullfile(files(k).folder, files(k).name)));
  if size(I, 3) == 3, I = mean(I, 3); end
  video(:, :, k) = I;
end
printf('%d frames of %dx%d px\n', size(video, 3), size(video, 2), size(video, 1));

%% PExFlow -------------------------------------------------------------------
[dX, dY] = PExFlow2D(video, window_size, poly_degree, pyr_levels, iterations, ...
                     sigma_int, frame_gap, 'Verbose', true);
F = pexflow_fields(dX, dY, dt * frame_gap, scale);

%% Save ----------------------------------------------------------------------
U = F.U; V = F.V; x = F.x; y = F.y; t = (0:size(U, 3) - 1) * dt;
convencion_V = 'positive = upwards (decreasing pixel row)';
save('-v7', 'pexflow_result.mat', 'U', 'V', 'x', 'y', 't', 'convencion_V');
printf('Saved pexflow_result.mat\n');

%% Quick look at the first field -------------------------------------------
figure('Color', 'w');
s = 16; [xs, ys] = meshgrid(1:s:size(U, 2), 1:s:size(U, 1));
imagesc(F.vorticity(:, :, 1)); axis image; colorbar; hold on;
quiver(xs, ys, U(1:s:end, 1:s:end, 1), -V(1:s:end, 1:s:end, 1), 'k');
title('vorticity [1/s] and velocity vectors (frame pair 1)');
