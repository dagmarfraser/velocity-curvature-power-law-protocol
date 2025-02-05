clear all
close all

addpath(genpath('functions'))
addpath(genpath('req'))

% Load appropriate data based on paramChoice
paramChoice = 2; % 1 for Maoz, 2 for Schaal
if paramChoice == 1
    load('Maoz_paperFigure.mat')
    plotTitle = {'Comparing Regression Methods', ...
                 'Constant Tangential Velocity + White Noise (Maoz et al., 2005)'};
else
    load('Schaal_paperFigure.mat')
    plotTitle = {'Comparing Filtering and Regression Methods', ...
                 'Constant Tangential Velocity + White Noise (Schaal & Sternad, 2001)'};
end

% Enhanced figure settings with larger text
set(groot, ...
    'defaultLineLineWidth', 2.0, ...
    'defaultAxesFontSize', 14, ...
    'defaultTextFontSize', 14, ...
    'defaultAxesLineWidth', 1.5, ...
    'defaultFigureColor', 'w', ...
    'defaultAxesBox', 'on', ...
    'defaultAxesTickDir', 'out');

% Define colors and line styles
colors = [    
    0.8500    0.3250    0.0980;  % Dark orange    
    0.0000    0.4470    0.7410;  % Blue
    0.4940    0.1840    0.5560   % Purple
];

lineStyles = {'-', '--', '-.'};  % Solid, dashed, dash-dot

% Create larger figure based on paramChoice
if paramChoice == 1
    % Create single plot for Maoz with larger dimensions
    fig = figure('Position', [50 50 1200 800]);
    
    % Get screen size and adjust if needed
    screen = get(0, 'ScreenSize');
    if screen(3) < 1200 || screen(4) < 800
        figPos = [(screen(3)-1200)/2 (screen(4)-800)/2 1200 800];
        set(fig, 'Position', figPos);
        set(fig, 'WindowStyle', 'normal');
    end
    
    % Use axes instead of tiledlayout for Maoz
    ax = axes;
    hold on
    
    % Store line handles for legend
    lineHandles = gobjects(1,3);
    
    % Plot data for regression types 3-5 only, flipping sign of beta values
    for regIdx = 3:5
        % Get data and calculate min/max envelope
        result = -reshape(resultDiff(:,regIdx), [], reps)';  % Flip sign here
        meanLine = nanmean(result);
        minLine = min(result, [], 1, 'omitnan');
        maxLine = max(result, [], 1, 'omitnan');
        
        % Create shaded min/max envelope
        fill([noise(1:noiseLength) fliplr(noise(1:noiseLength))], ...
             [maxLine fliplr(minLine)], ...
             colors(regIdx-2,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
         
        % Plot mean line with different line styles
        lineHandles(regIdx-2) = plot(noise(1:noiseLength), meanLine, ...
             lineStyles{regIdx-2}, 'Color', colors(regIdx-2,:), 'LineWidth', 2.5);
    end
    
else
    % Create double plot for Schaal with larger dimensions
    fig = figure('Position', [50 50 1800 900]);
    
    % Get screen size and adjust if needed
    screen = get(0, 'ScreenSize');
    if screen(3) < 1800 || screen(4) < 900
        figPos = [(screen(3)-1800)/2 (screen(4)-900)/2 1800 900];
        set(fig, 'Position', figPos);
        set(fig, 'WindowStyle', 'normal');
    end
    
    t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    % Store line handles for each subplot
    lineHandles1 = gobjects(1,3);
    lineHandles2 = gobjects(1,3);
    
    % First plot - Butterworth
    ax1 = nexttile;
    hold on
    for regIdx = 3:5
        result = -reshape(resultDiffBW(:,regIdx), [], reps)';  % Flip sign here
        meanLine = nanmean(result);
        minLine = min(result, [], 1, 'omitnan');
        maxLine = max(result, [], 1, 'omitnan');
        
        fill([noise(1:noiseLength) fliplr(noise(1:noiseLength))], ...
             [maxLine fliplr(minLine)], ...
             colors(regIdx-2,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        lineHandles1(regIdx-2) = plot(noise(1:noiseLength), meanLine, ...
             lineStyles{regIdx-2}, 'Color', colors(regIdx-2,:), 'LineWidth', 2.5);
    end
    title({'Butterworth Filter', 'Exhibiting Spurious 1/3 Power Law'}, 'FontSize', 16)
    
    % Second plot - Savitzky-Golay
    ax2 = nexttile;
    hold on
    for regIdx = 3:5
        result = -reshape(resultSG(:,regIdx), [], reps)';  % Flip sign here
        meanLine = nanmean(result);
        minLine = min(result, [], 1, 'omitnan');
        maxLine = max(result, [], 1, 'omitnan');
        
        fill([noise(1:noiseLength) fliplr(noise(1:noiseLength))], ...
             [maxLine fliplr(minLine)], ...
             colors(regIdx-2,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        lineHandles2(regIdx-2) = plot(noise(1:noiseLength), meanLine, ...
             lineStyles{regIdx-2}, 'Color', colors(regIdx-2,:), 'LineWidth', 2.5);
    end
    title({'Savitzky-Golay Filter', 'Avoiding Spurious 1/3 Power Law'}, 'FontSize', 16)
end

% Add common elements to all axes
ax = findall(gcf, 'type', 'axes');
for i = 1:length(ax)
    axes(ax(i));
    
    % Calculate appropriate y-axis limits
    if i == 1
        allLines = findobj(ax(i), 'Type', 'line');
        yData = cell2mat(get(allLines, 'YData'));
        yMax = max(yData(:));
        yMin = min(yData(:));
    end
    yRange = yMax - yMin;
    extraSpace = yRange * 0.15;
    upperBound = 1/3 * 1.1;  % Calculate bounds first to use in limits
    lowerBound = 1/3 * 0.9;
    
    plotYMin = min([yMin - extraSpace, 0, lowerBound - 0.05]);  % Ensure zero and lower bound are visible
    plotYMax = max([yMax + extraSpace, upperBound + 0.05]);  % Ensure upper bound is visible with padding
    
    % Add reference lines with enhanced visibility
    yline(1/3, '-', 'LineWidth', 2.5, 'Color', [0.2 0.2 0.8]);  % Changed from -1/3
    text(min(noise), 1/3, '  β = 1/3', 'VerticalAlignment', 'bottom', ...  % Changed from -1/3
         'FontWeight', 'bold', 'Color', [0.2 0.2 0.8], 'FontSize', 14, 'Interpreter', 'tex');
    
    % Add ±10% reference lines
    upperBound = 1/3 * 1.1;  % Changed bounds
    lowerBound = 1/3 * 0.9;
    
    plot([min(noise) max(noise)], [upperBound upperBound], ':', ...
         'Color', [0.2 0.2 0.8], 'LineWidth', 2);
    plot([min(noise) max(noise)], [lowerBound lowerBound], ':', ...
         'Color', [0.2 0.2 0.8], 'LineWidth', 2);
    
    text(min(noise), upperBound, '  +10%', 'VerticalAlignment', 'bottom', ...
         'Color', [0.2 0.2 0.8], 'FontSize', 14, 'Interpreter', 'tex');
    text(min(noise), lowerBound, '  -10%', 'VerticalAlignment', 'top', ...
         'Color', [0.2 0.2 0.8], 'FontSize', 14, 'Interpreter', 'tex');
    
    % Customize appearance
    grid on
    box on
    xlabel('Noise Magnitude (mm)', 'FontSize', 16)
    ylabel('Power Law β Value', 'FontSize', 16, 'Interpreter', 'tex')  % Updated label
    
    % Set y-axis limits
    ylim([plotYMin plotYMax])
    xlim([min(noise)-0.5 max(noise)])
    
    % Increase tick label sizes
    set(gca, 'FontSize', 14)
end

% Add main title
if paramChoice == 1
    title(plotTitle, 'FontSize', 18);
else
    if verLessThan('matlab', '9.11')
        sgtitle(plotTitle, 'FontSize', 18);
    else
        title(t, plotTitle, 'FontSize', 18);
    end
end

% Add legends and annotations with improved spacing
legendStrings = {'Linear Regression', ...
                'Levenberg-Marquardt', ...
                'Iteratively Reweighted'};

if paramChoice == 1
    leg = legend(lineHandles, legendStrings, ...
                'Box', 'on', ...
                'EdgeColor', [0.7 0.7 0.7], ...
                'FontSize', 14, ...
                'Location', 'eastoutside');
    
    annotation('textbox', [0.15 0.85 0.3 0.1], ...
        'String', 'Shaded areas represent min/max envelope', ...
        'EdgeColor', 'none', 'FitBoxToText', 'on', 'FontSize', 14);
else
    leg1 = legend(ax1, lineHandles1, legendStrings, ...
                 'Box', 'on', ...
                 'EdgeColor', [0.7 0.7 0.7], ...
                 'FontSize', 14);
    leg1Pos = get(ax1, 'Position');
    set(leg1, 'Position', [leg1Pos(1) + leg1Pos(3)*0.4, ...
                          leg1Pos(2) + leg1Pos(4)*0.4, ...
                          0.2, 0.2]);
    
    leg2 = legend(ax2, lineHandles2, legendStrings, ...
                 'Box', 'on', ...
                 'EdgeColor', [0.7 0.7 0.7], ...
                 'FontSize', 14);
    leg2Pos = get(ax2, 'Position');
    set(leg2, 'Position', [leg2Pos(1) + leg2Pos(3)*0.4, ...
                          leg2Pos(2) + leg2Pos(4)*0.4, ...
                          0.2, 0.2]);
    
    annotation('textbox', [0.2 0.85 0.25 0.1], ...
        'String', 'Shaded areas represent min/max envelope', ...
        'EdgeColor', 'none', 'FitBoxToText', 'on', 'FontSize', 14);
    
    annotation('textbox', [0.65 0.85 0.25 0.1], ...
        'String', 'Shaded areas represent min/max envelope', ...
        'EdgeColor', 'none', 'FitBoxToText', 'on', 'FontSize', 14);
end

% Save figure if needed
if exist('saveImage', 'var') && saveImage
    folderSave = pwd;
    folderSave = [folderSave(1:end-3) 'figures'];
    
    fig = gcf;
    plotedit(fig, 'on');
    disp('PRESS ANY KEY ONCE PLOT ELEMENTS MOVED TO FINAL POSITION!')
    pause
    plotedit(fig, 'off');
    
    if paramChoice == 1
        figFilename = 'Maoz_revised';
    else
        figFilename = 'Schaal_revised';
    end
    
    % Save with higher resolution
    exportgraphics(fig, fullfile(folderSave, [figFilename, '.png']), ...
        'Resolution', 1200, ...
        'BackgroundColor', 'white', ...
        'ContentType', 'vector');
        
    exportgraphics(fig, fullfile(folderSave, [figFilename, '.eps']), ...
        'Resolution', 1200);
end