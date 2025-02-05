% Enhanced visualization combining beta analysis and shape visualization
% with publication-quality settings

clear all
close all

addpath(genpath('functions'))
addpath(genpath('req'))

saveImage = 1;
% Publication-quality figure settings
set(groot, ...
    'defaultLineLineWidth', 2.0, ...
    'defaultAxesFontSize', 14, ...
    'defaultTextFontSize', 14, ...
    'defaultAxesLineWidth', 1.5, ...
    'defaultFigureColor', 'w', ...
    'defaultAxesBox', 'on', ...
    'defaultAxesTickDir', 'out');

% Create figure with appropriate size for both plots
fig = figure('Position', [50 50 1800 900]);

% Get screen size and adjust if needed
screen = get(0, 'ScreenSize');
if screen(3) < 1800 || screen(4) < 900
    figPos = [(screen(3)-1800)/2 (screen(4)-900)/2 1800 900];
    set(fig, 'Position', figPos);
    set(fig, 'WindowStyle', 'normal');
end

% Angular frequencies (σ)
sigmaList = [0.01 2/5 3/5 2/3 4/5 4/3 3/2 2 5/2 3 4 6];
sigmaListCell = {'0.01' '2/5' '3/5' '2/3' '4/5' '4/3' '3/2' '2' '5/2' '3' '4' '6'};

% Left subplot for beta analysis - adjusted position
pos_left = [0.08 0.1 0.38 0.8];  % Moved right edge closer to center
subplot('Position', pos_left)

% Dense list for continuous line
sigma_dense = logspace(-2, log10(6), 1000);

% Beta calculation function based on angular frequency (σ)
beta = @(sigma) (2/3) * ((1 + sigma.^2/2) ./ (1 + sigma.^2 + sigma.^4/15));

beta_dense = beta(sigma_dense);
beta_points = beta(sigmaList);

hold on
% Plot continuous line with enhanced visibility and store handle
h1 = plot(sigma_dense, beta_dense, 'LineWidth', 2.5, 'Color', [0.2 0.6 0.8], 'DisplayName', '\beta(\sigma) curve');
% Plot specific points with larger markers and store handle
h2 = scatter(sigmaList, beta_points, 100, 'ro', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', 'Points of interest');

% Reference line
yline(1/3, '--k', '\beta = 1/3', 'LineWidth', 1.5, 'FontSize', 14, 'LabelHorizontalAlignment', 'left');

% Enhanced axes properties
grid on
box on
set(gca, 'XScale', 'log')
xlabel('Angular Frequency (\sigma)', 'FontSize', 16, 'Interpreter', 'tex');
ylabel('Power Law Exponent (\beta)', 'FontSize', 16, 'Interpreter', 'tex');
title({'\newline', 'Power Law Exponent vs Angular Frequency', ''}, 'FontSize', 18);

% Add point labels with improved positioning
for i = 1:length(sigmaList)
    text(sigmaList(i), beta_points(i) + 0.02, sigmaListCell{i}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 12);
end

% Refined axis limits
set(gca, 'XLim', [0.008 7], 'YLim', [0.1 0.9]);

% Add legend with explicit strings
legend([h1, h2], {'\beta(\sigma) curve', 'Points of interest'}, 'Location', 'northeast', 'FontSize', 14);

% Formula display with enhanced visibility
text(0.015, 0.8, {'$\mathbf{\beta(\sigma) = \frac{2}{3} \cdot \frac{1 + \frac{\sigma^2}{2}}{1 + \sigma^2 + \frac{\sigma^4}{15}}}$'}, ...
    'FontSize', 28, 'Interpreter', 'latex', 'HorizontalAlignment', 'left');

% Custom grid improvements
set(gca, 'XGrid', 'off', 'YGrid', 'on')
grid minor

% Add custom major grid lines
axisLims = axis;
for sigma = [0.01 0.1 1 10]
    if sigma >= axisLims(1) && sigma <= axisLims(2)
        line([sigma sigma], [axisLims(3) axisLims(4)], ...
            'Color', [0.8 0.8 0.8], 'LineStyle', '-', 'LineWidth', 0.5);
    end
end

% Right side subplots for shapes - adjusted position
pos_right = [0.52 0.1 0.38 0.8];  % Moved left edge closer to center
margin = 0.015;  % Reduced margin between subplots
nRows = 3;
nCols = 4;
width = (pos_right(3) - (nCols+1)*margin) / nCols;
height = (pos_right(4) - (nRows+1)*margin) / nRows;

% Shape generation parameters
epsilon = 1.3;
theta0 = 0;

% Generate and plot shapes
for i = 1:length(sigmaList)
    % Calculate row and column indices
    row = ceil(i/nCols);
    col = mod(i-1, nCols) + 1;
    
    % Calculate position for this subplot
    left = pos_right(1) + margin + (col-1)*(width + margin);
    bottom = pos_right(2) + pos_right(4) - (row*height + margin*row);
    pos = [left bottom width height];
    
    % Create subplot and plot shape
    ax = subplot('Position', pos);
    [x, y] = pureCurveGenerator(0, 0, epsilon, sigmaList(i), theta0, 10, 200);
    plot(x, y, 'LineWidth', 2, 'Color', [0.2 0.6 0.8]);
    axis equal off
    title(sprintf('\\sigma = %s', sigmaListCell{i}), ...
        'FontSize', 14, 'Interpreter', 'tex');
end

% Add title for the shapes section with blank line
axes('Position', [pos_right(1) pos_right(2)+pos_right(4) pos_right(3) 0.02], 'Visible', 'off');
text(0.5, 0, {'\newline', 'Elementary Shapes of Angular Frequency \sigma'}, ...
    'FontSize', 18, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom');

% Overall figure title
sgtitle('Relationship Between Angular Frequency, Power Law Exponent, and Shape', ...
    'FontSize', 20, 'FontWeight', 'bold');

% Save figure if needed
if exist('saveImage', 'var') && saveImage
    folderSave = pwd;
    folderSave = [folderSave(1:end-3) 'figures'];
    
    plotedit(fig, 'on');
    disp('PRESS ANY KEY ONCE PLOT ELEMENTS MOVED TO FINAL POSITION - PROPERTY INSPECTOR FOR LEGEND!!')
    pause
    plotedit(fig, 'off');
    
    % Save with high resolution
    exportgraphics(fig, fullfile(folderSave, 'beta_shapes_combined.png'), ...
        'Resolution', 1200, ...
        'BackgroundColor', 'white', ...
        'ContentType', 'vector');
    
    exportgraphics(fig, fullfile(folderSave, 'beta_shapes_combined.eps'), ...
        'Resolution', 1200);
end