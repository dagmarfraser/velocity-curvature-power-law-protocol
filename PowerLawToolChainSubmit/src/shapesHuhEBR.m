% consider https://github.com/adam-matic/KinematicCognition/blob/master/app/src/main/java/com/example/kinematiccognition/PureCurveGenerator.kt

% add a graph of the power law against beta
clear all
close all

addpath(genpath('functions'))
addpath(genpath('req'))

% nu = 2 %frequency
epsilon = 1 % amplitude
theta0 = 0 % rotation
num_points = 1000

epsilons = 1.3 ;linspace(0, 2, 100);
nuList = [0.01 2/5 3/5 2/3 4/5 4/3 3/2 2 5/2 3 4 6] ;
nuListCell = {'0.01' '2/5' '3/5' '2/3' '4/5' '4/3' '3/2' '2' '5/2' '3' '4' '6'} ;
fig = figure
t = tiledlayout(3,4);
t.TileSpacing = 'compact';

for nus = 1:length(nuList)
    for its = 1:length(epsilons)

        [x,y] = pureCurveGenerator(0, 0, epsilons(its), nuList(nus), theta0, 10, 200);
        nexttile
        plot(x, y, 'LineWidth',3);
        axis equal;
        axis off;
        title(sprintf('\\nu = %s', nuListCell{nus}), 'FontSize', 14, 'Interpreter', 'tex');
      %  title(sprintf('\\nu = %s', nuListCell{nus}), 'FontSize', 14, 'Interpreter', 'latex');
      %  title(sprintf('\nu = %s', nuListCell{nus}), 'FontSize', 28);
%        title(sprintf('Elementary Shape: ν = %s, ε = %.2f, θ_0 = %.2f', nuListCell{nus}, epsilons(its), theta0));
        drawnow
    end
end

% Adjust layout
title(t, 'Elementary Shapes of Angular Frequency \nu', 'FontSize', 16, 'FontWeight', 'bold');
disp('press any key to save')
pause
% Save as high-quality EPS
print(fig, 'elementary_shapes', '-depsc', '-r600', '-painters');
