% Enhanced visualization combining beta analysis and shape visualization
% with publication-quality settings

clear all
close all

addpath(genpath('functions'))
addpath(genpath('req'))

figFilename11 = 'beta_shapes_combined';
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

% Angular frequencies (φ)
phiList = [0.01 2/5 3/5 2/3 4/5 4/3 3/2 2 5/2 3 4 6];
phiListCell = {'0.01' '2/5' '3/5' '2/3' '4/5' '4/3' '3/2' '2' '5/2' '3' '4' '6'};

% Left subplot for beta analysis - adjusted position
pos_left = [0.08 0.1 0.38 0.8];  % Moved right edge closer to center
ax_left = subplot('Position', pos_left);

% Dense list for continuous line
phi_dense = logspace(-2, log10(6), 1000);

% Beta calculation function based on angular frequency (φ)
beta = @(phi) (2/3) * ((1 + phi.^2/2) ./ (1 + phi.^2 + phi.^4/15));

beta_dense = beta(phi_dense);
beta_points = beta(phiList);

hold on
% Plot continuous line with enhanced visibility and store handle
h1 = plot(phi_dense, beta_dense, 'LineWidth', 2.5, 'Color', [0.2 0.6 0.8], 'DisplayName', '\beta(\phi) curve');
% Plot specific points with larger markers and store handle
h2 = scatter(phiList, beta_points, 80, 'ro', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', 'Points of interest');

% Reference line
yline(1/3, '--k', '\beta = 1/3', 'LineWidth', 2, 'FontSize', 16, 'LabelHorizontalAlignment', 'left');

% Enhanced axes properties
grid on
box on
set(ax_left, 'XScale', 'log')
xlabel('Angular Frequency (\phi)', 'FontSize', 18, 'Interpreter', 'tex');
ylabel('Power Law Exponent (\beta)', 'FontSize', 18, 'Interpreter', 'tex');
title('\beta(\phi) vs \phi', 'FontSize', 20);

% Add point labels with improved positioning
for i = 1:length(phiList)
    text(phiList(i), beta_points(i) + 0.02, phiListCell{i}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 14);
end

% Refined axis limits
xlim([0.008 7]);
ylim([0.1 0.9]);

% Add legend with explicit strings
legend([h1, h2], {'\beta(\phi) curve', 'Points of interest'}, 'Location', 'northeast', 'FontSize', 16);

% Custom grid improvements
set(ax_left, 'XGrid', 'off', 'YGrid', 'on')
grid minor

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
for i = 1:length(phiList)
    % Calculate row and column indices
    row = ceil(i/nCols);
    col = mod(i-1, nCols) + 1;
    
    % Calculate position for this subplot
    left = pos_right(1) + margin + (col-1)*(width + margin);
    bottom = pos_right(2) + pos_right(4) - (row*height + margin*row);
    pos = [left bottom width height];
    
    % Create subplot and plot shape
    ax = subplot('Position', pos);
    [x, y] = pureCurveGenerator(0, 0, epsilon, phiList(i), theta0, 10, 200);
    plot(x, y, 'LineWidth', 2, 'Color', [0.2 0.6 0.8]);
    axis equal off
    title(sprintf('\\phi = %s', phiListCell{i}), ...
        'FontSize', 16, 'Interpreter', 'tex');
    
    % Apply some styling to shape subplots
    set(ax, 'FontName', 'Arial', 'FontSize', 14);
end

% Add title for the shapes section with blank line
axes('Position', [pos_right(1) pos_right(2)+pos_right(4) pos_right(3) 0.02], 'Visible', 'off');
text(0.5, 0, {'', 'Elementary Shapes of Angular Frequency \phi'}, ...
    'FontSize', 20, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom');

% Overall figure title
sgtitle('Relationship Between Angular Frequency, Power Law Exponent, and Shape', ...
    'FontSize', 24, 'FontWeight', 'bold');

% Export with standardized settings
if saveImage
    % Enable manual adjustment before export
    exportPublicationFigure(fig, figFilename11, ...
        'Width', 24, 'Height', 18, ...
        'Formats', {'png', 'eps'}, ...
        'Resolution', 1200, ...
        'ManualAdjust', true);  % Add manual adjustment parameter
end