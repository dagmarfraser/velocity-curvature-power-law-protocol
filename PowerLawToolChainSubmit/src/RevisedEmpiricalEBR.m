clear all
close all
load

figFilename11 = 'Zarandi_Empirical_TOST';

% Process beta and VGF data - flip sign of beta values
for diffType = 2:3
    for regressType = 3:5
        betaArray = -beta(:,:,diffType,regressType); % Flip sign here
        betaSerial = betaArray(:);
        betaSerialCell{diffType, regressType} = betaSerial;
        betaMean(diffType, regressType) = nanmean(betaSerial);
        betaStd(diffType, regressType) = nanstd(betaSerial);
    end
end

% Create grouping variables
g1 = ones(size(betaSerialCell{2,3})) * 1;
g2 = ones(size(betaSerialCell{2,4})) * 2;
g3 = ones(size(betaSerialCell{2,5})) * 3;
g4 = ones(size(betaSerialCell{3,3})) * 4;
g5 = ones(size(betaSerialCell{3,4})) * 5;
g6 = ones(size(betaSerialCell{3,5})) * 6;

% Test for equivalence
target = 1/3;  % Changed from -1/3 to 1/3
margin = 0.1;
methodNames = {'BW LR', 'BW LMLS', 'BW IRLS', 'SG LR', 'SG LMLS', 'SG IRLS'};
equivalence = zeros(6,1);

allData = {betaSerialCell{2,3}, betaSerialCell{2,4}, betaSerialCell{2,5}, ...
           betaSerialCell{3,3}, betaSerialCell{3,4}, betaSerialCell{3,5}};

lower_bound = target * (1 - margin);  % Changed order of bounds
upper_bound = target * (1 + margin);

% Test each method
for i = 1:6
    data = allData{i};
    data = data(~isnan(data));
    sample_mean = mean(data);
    equivalence(i) = (sample_mean >= lower_bound) && (sample_mean <= upper_bound);
    
    % Display diagnostic information
    disp(['Method ' num2str(i) ' (' methodNames{i} '):']);
    disp(['  Mean: ' num2str(sample_mean)]);
    disp(['  Lower bound: ' num2str(lower_bound)]);
    disp(['  Upper bound: ' num2str(upper_bound)]);
    disp(['  Equivalent: ' num2str(equivalence(i))]);
end

% Create larger figure
fig = figure(11);

% Set larger dimensions
screen = get(0, 'ScreenSize');
figWidth = 1800;
figHeight = 1200;

% Position figure with screen size check
if screen(3) < figWidth || screen(4) < figHeight
    figPos = [(screen(3)-figWidth)/2 (screen(4)-figHeight)/2 figWidth figHeight];
    set(fig, 'Position', figPos);
    set(fig, 'WindowStyle', 'normal');
else
    set(fig, 'Position', [50 50 figWidth figHeight]);
end

% Create box plot with improved styling
h = boxplot([betaSerialCell{2,3}; betaSerialCell{2,4}; betaSerialCell{2,5}; ...
    betaSerialCell{3,3}; betaSerialCell{3,4}; betaSerialCell{3,5}],  ...
    [g1;g2;g3;g4;g5;g6], ...
    'Notch', 'on', ...
    'Labels', methodNames, ...
    'Colors', [0.2 0.2 0.7], ...
    'Width', 0.7);

% Enhanced boxplot appearance
set(findobj(gca,'type','line'), 'LineWidth', 2)
set(findobj(gca,'tag','Median'), 'LineWidth', 2.5, 'Color', [0.8 0 0])
set(findobj(gca,'tag','Outliers'), 'MarkerEdgeColor', [0.2 0.2 0.7], ...
    'MarkerFaceColor', [0.8 0.8 1], 'MarkerSize', 10)

% Calculate appropriate y-axis limits that include the bounds
yMax = max(cellfun(@max, allData));
yMin = min(cellfun(@min, allData));
yRange = yMax - yMin;
extraSpace = yRange * 0.15;  % 15% extra space for text and markers
plotYMin = min(yMin - extraSpace, lower_bound - extraSpace/2);  % Ensure lower bound is visible
plotYMax = max(yMax + extraSpace/2, upper_bound + extraSpace/2);  % Ensure upper bound is visible

% Add grid and reference lines
grid on
yl = yline([1/3], '-.r', 'β = 1/3', 'LineWidth', 2.5, 'FontSize', 16, 'FontWeight', 'bold');  % Changed from -1/3
yl.LabelHorizontalAlignment = 'left';
yl.LabelVerticalAlignment = 'bottom';

% Add equivalence bounds
yl_lower = yline([lower_bound], ':k', ['Lower bound: ' num2str(lower_bound, '%.3f')], 'LineWidth', 1.5, 'FontSize', 14);
yl_upper = yline([upper_bound], ':k', ['Upper bound: ' num2str(upper_bound, '%.3f')], 'LineWidth', 1.5, 'FontSize', 14);
yl_lower.LabelHorizontalAlignment = 'left';
yl_upper.LabelHorizontalAlignment = 'left';

% Position significance markers
grid off
starPos = plotYMin + extraSpace/2;  % Position markers above the bottom of the plot

% Create background boxes for significance markers
for i = 1:6
    rectangle('Position', [i-0.25, starPos-0.02*yRange, 0.5, 0.04*yRange], ...
              'FaceColor', 'white', ...
              'EdgeColor', 'none');
end

% Add significance markers
for i = 1:6
    if equivalence(i)
        text(i, starPos, 'equiv', 'HorizontalAlignment', 'center', ...
             'FontSize', 18, 'FontWeight', 'bold', 'Color', [0 0.5 0]);
    else
        text(i, starPos, 'non-equiv', 'HorizontalAlignment', 'center', ...
             'FontSize', 18, 'FontWeight', 'bold', 'Color', [0.8 0 0]);
    end
end

% Restore grid and adjust layer order
grid on
set(gca, 'Layer', 'bottom')

% Set y-axis limits to include bounds and markers
ylim([plotYMin plotYMax])

% Enhanced titles and labels
title({'Comparison of Power Law β Values Across Analysis Methods', ...
       'Using Empirical Ellipse Tracing Data (Zarandi et al. 2023)'}, ...
    'FontSize', 20, 'FontWeight', 'bold')

ylabel('Power Law β Value', 'FontSize', 18, 'FontWeight', 'bold')
xlabel({' ', ...
    'Filter Types: BW (Butterworth), SG (Savitzky-Golay)', ...
    'Regression Types: LR (Linear), LMLS (Levenberg-Marquardt), IRLS (Iteratively Reweighted)', ...
    ' '}, 'FontSize', 16)

% Update box plot explanation
dim = [.15 .7 .2 .2];
annotation('textbox', dim, 'String', ...
    {'Box Plot Elements:', ...
     '• Box: 25th to 75th percentiles', ...
     '• Center line: median', ...
     '• Notches: 95% CI of median', ...
     '• Whiskers: data range (excluding outliers)', ...
     '• Points: outliers'}, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', ...
    'FontSize', 14);

% Add equivalence testing explanation
dim = [.7 .7 .2 .15];
annotation('textbox', dim, 'String', ...
    {'Equivalence Testing:', ...
     'equiv: mean within ±10% of 1/3', ...  % Changed from -1/3
     'non-equiv: mean outside ±10% bounds'}, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', ...
    'FontSize', 14);

% Enhance axes properties
ax = gca;
ax.FontSize = 14;
ax.LineWidth = 1.5;
ax.Box = 'on';

% Export the figure if requested
if exist('saveImage', 'var') && saveImage
    folderSave = pwd;
    folderSave = [folderSave(1:end-3) 'figures'];
    
    % Enable manual adjustment mode
    plotedit(fig, 'on');
    disp('PRESS ANY KEY ONCE PLOT ELEMENTS MOVED TO FINAL POSITION!')
    pause
    plotedit(fig, 'off');
    
    % Save with high resolution
    exportgraphics(fig, fullfile(folderSave, [figFilename11, '.png']), ...
        'Resolution', 1200, ...
        'BackgroundColor', 'white', ...
        'ContentType', 'vector');
    
    exportgraphics(fig, fullfile(folderSave, [figFilename11, '.eps']), ...
        'Resolution', 1200);
end