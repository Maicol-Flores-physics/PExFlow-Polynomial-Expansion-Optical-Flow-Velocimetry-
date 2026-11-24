function F = pexflow_fields(dX, dY, dt, scale)
%PEXFLOW_FIELDS  Physical velocity and derived fields from PExFlow2D output.
%
%   F = pexflow_fields(dX, dY, dt, scale)
%
%   dX, dY  M x N x T displacements [px] returned by PExFlow2D
%   dt      time between the two frames of each pair [s] (dt * frame_gap)
%   scale   calibration [m/px] (or any length unit per pixel)
%
%   F is a struct with fields (all M x N x T, physical convention: x to
%   the right, y UPWARDS):
%     U, V        velocity components            [scale/s]
%     speed       |u|                            [scale/s]
%     vorticity   dV/dx - dU/dy                  [1/s]
%     divergence  dU/dx + dV/dy                  [1/s]
%     angle       atan2(V, U)                    [deg]
%     x, y        coordinates of the columns / rows [scale]
%
%   MATLAB R2018a+ / GNU Octave >= 6.
if nargin < 3 || isempty(dt),    dt = 1;    end
if nargin < 4 || isempty(scale), scale = 1; end
[M, N, T] = size(dX);
F.U = double(dX) * (scale / dt);
F.V = -double(dY) * (scale / dt);          % image rows grow downwards
F.speed = hypot(F.U, F.V);
F.vorticity  = zeros(M, N, T);
F.divergence = zeros(M, N, T);
for t = 1:T
    [dUdx, dUdy_img] = gradient(F.U(:, :, t), scale);
    [dVdx, dVdy_img] = gradient(F.V(:, :, t), scale);
    F.vorticity(:, :, t)  = dVdx + dUdy_img;   % dU/dy = -dU/drow
    F.divergence(:, :, t) = dUdx - dVdy_img;
end
F.angle = atan2(F.V, F.U) * 180 / pi;
F.x = (0:N - 1) * scale;
F.y = (0:M - 1) * scale;
end
