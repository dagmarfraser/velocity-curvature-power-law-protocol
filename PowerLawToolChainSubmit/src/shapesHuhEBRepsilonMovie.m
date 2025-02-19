% consider https://github.com/adam-matic/KinematicCognition/blob/master/app/src/main/java/com/example/kinematiccognition/PureCurveGenerator.kt

% add a graph of the power law against beta
clear all
close all

addpath(genpath('functions'))
addpath(genpath('req'))

% nu = 2 %frequency
theta0 = 0 % rotation
num_points = 1000

nuList = [0.01 2/5 3/5 2/3 4/5 4/3 3/2 2 5/2 3 4 6];
nuListCell = {'0.01' '2/5' '3/5' '2/3' '4/5' '4/3' '3/2' '2' '5/2' '3' '4' '6'};

% Create animation parameters
epsilon_min = 0;
epsilon_max = 2;
n_frames = 60;
epsilons = linspace(epsilon_min, epsilon_max, n_frames);

% Setup figure
fig = figure('Position', [100, 100, 1200, 800]);
t = tiledlayout(3, 4);
t.TileSpacing = 'compact';
title(t, 'Elementary Shapes of Angular Frequency \nu', 'FontSize', 16, 'FontWeight', 'bold');

% Create VideoWriter object
v = VideoWriter('epsilon_animation.mp4', 'MPEG-4');
v.FrameRate = 10;
v.Quality = 100;
open(v);

% Animation loop
for frame = 1:n_frames
    epsilon = epsilons(frame);
    
    % Clear previous plots but keep the figure
    clf;
    t = tiledlayout(3, 4);
    t.TileSpacing = 'compact';
    
    for nus = 1:length(nuList)
        [x, y] = pureCurveGenerator(0, 0, epsilon, nuList(nus), theta0, 10, 200);
        nexttile
        plot(x, y, 'LineWidth', 3);
        axis equal;
        axis off;
        title(sprintf('\\nu = %s', nuListCell{nus}), 'FontSize', 14, 'Interpreter', 'tex');
    end
    
    % Add main title with current epsilon value
    title(t, sprintf('Elementary Shapes: \\nu varying, \\epsilon = %.2f', epsilon), 'FontSize', 16, 'FontWeight', 'bold', 'Interpreter', 'tex');
    
    % Capture the frame
    drawnow;
    frame_data = getframe(fig);
    writeVideo(v, frame_data);
end

% Close the video file
close(v);
disp('Animation saved as epsilon_animation.mp4');