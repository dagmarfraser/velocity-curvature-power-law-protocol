% PlotCitesR2R_Updated.m
% Updated version using standardized styling functions for publication

clear all
close all

addpath(genpath('functions'))
addpath(genpath('req'))

% Load the data
load('LacquanitiCites_R2.mat')

% Get consistent styling parameters
[colors, lineStyles] = getConsistentStyles();

% Figure 1: Estimated Citations per Year
fig1 = figure(1);
plot(VarName3, VarName1, 'Color', colors(1,:), 'LineStyle', lineStyles{1}, ...
    'Marker', '.', ...
    'MarkerFaceColor', 'white');
ylabel('Citations from Google Scholar');
xlabel('Year');
title('Citations of Lacquaniti et al. (1983) per year');
grid off;

% Set x-axis limit from 1980 to 2025
xlim([1980 2025]);

% Create more consistent tick spacing with major and minor ticks
% Major ticks every 5 years
majorTicks = 1980:5:2025;
set(gca, 'XTick', majorTicks);

% Add minor ticks for individual years
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 1980:1:2025;

% Format tick labels and rotate for better spacing
xticklabels(cellstr(num2str(majorTicks', '%d')));
set(gca, 'XTickLabelRotation', 45);  % Rotate labels by 45 degrees

% % Add a vertical line at 1983 to highlight the publication year
% hold on;
% xline(1983, '--k', 'Publication Year', 'LineWidth', 1.5, 'Alpha', 0.7, 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');
% hold off;

% Apply publication styling
setPublicationStyle(fig1);

% Export with standardized settings
exportPublicationFigure(fig1, 'citations_per_year', ...
    'Width', 16, 'Height', 12, ...
    'Formats', {'png', 'eps'});

% Figure 2: Estimated Total Citations
fig2 = figure(2);
plot(flip(VarName3), cumsum(flip(VarName1)), 'Color', colors(1,:), 'LineStyle', lineStyles{1}, ...
    'Marker', '.', ...
    'MarkerFaceColor', 'white');
ylabel('Est. Total Citations');
xlabel('Year');
title('Total Citations of Lacquaniti et al. (1983)');

% Set x-axis limit from 1980 to 2025
xlim([1980 2025]);

% Create more consistent tick spacing with major and minor ticks
% Major ticks every 5 years
majorTicks = 1980:5:2025;
set(gca, 'XTick', majorTicks);

% Add minor ticks for individual years
ax = gca;
ax.XAxis.MinorTick = 'on';
ax.XAxis.MinorTickValues = 1980:1:2025;

% Format tick labels and rotate for better spacing
xticklabels(cellstr(num2str(majorTicks', '%d')));
set(gca, 'XTickLabelRotation', 45);  % Rotate labels by 45 degrees

% % Add a vertical line at 1983 to highlight the publication year
% hold on;
% xline(1983, '--k', 'Publication Year', 'LineWidth', 1.5, 'Alpha', 0.7, 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');
% hold off;

% Apply publication styling
setPublicationStyle(fig2);

% Export with standardized settings
exportPublicationFigure(fig2, 'total_citations', ...
    'Width', 16, 'Height', 12, ...
    'Formats', {'png', 'eps'});