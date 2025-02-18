% code for Exp Brain Res review paper of Fraser et al., 2024
% As outlined in Huh 2015, and in Matic and Gomez-Marin 2019
% https://github.com/adam-matic/KinematicCognition/blob/master/app/src/main/java/com/example/kinematiccognition/PureCurveGenerator.kt
% adapted from pyhton to MATLAB
%
% Created July 2024
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk

function [xs, ys] = PureCurveGenerator(x0, y0, eps, v, th0, N, scale)
    % Set default values if not provided
    if nargin < 1, x0 = 0; end
    if nargin < 2, y0 = 0; end
    if nargin < 3, eps = 1.2; end
    if nargin < 4, v = 2; end
    if nargin < 5, th0 = 0; end
    if nargin < 6, N = 5; end
    if nargin < 7, scale = 200; end

    % Initialize variables
    xs = [];
    ys = [];
    th = th0;
    dtheta = 0.005;
    x = x0;
    y = y0;

    % Main loop
    while th < (th0 + pi * N)

        dx = dtheta * cos(th) * exp(eps * sin(v * (th - th0)));
        x = x + scale * dx;

        dy = dtheta * sin(th) * exp(eps * sin(v * (th - th0)));
        y = y + scale * dy;

        th = th + dtheta;

        xs(end+1) = x;
        ys(end+1) = y;

    end
end