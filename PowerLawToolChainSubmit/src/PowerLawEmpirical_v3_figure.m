% RevisedEmpiricalEBR_Updated.m
% Updated version using standardized styling functions for publication

clear all
close all
% Add required paths
addpath(genpath('functions'));
addpath(genpath('req'));

load % Load default matlab.mat with beta data

figFilename11 = 'Zarandi_Empirical_TOST';

% Get consistent styling parameters
[colors, lineStyles] = getConsistentStyles();

% Process beta and VGF data - flip sign of beta valuestextbox
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

% Create figure
fig = figure();

% Create box plot
h = boxplot([betaSerialCell{2,3}; betaSerialCell{2,4}; betaSerialCell{2,5}; ...
    betaSerialCell{3,3}; betaSerialCell{3,4}; betaSerialCell{3,5}],  ...
    [g1;g2;g3;g4;g5;g6], ...
    'Notch', 'on', ...
    'Labels', methodNames, ...
    'Colors', colors(1,:), ...
    'Width', 0.7);

title(gca, {'Comparison of Power Law β Values Across Analysis Methods', ...
            'Using Empirical Ellipse Tracing Data (Zarandi et al. 2023)'}, ...
            'FontSize', 16, ...
            'FontWeight', 'bold');

% Increase box plot line widths for better visibility in print
boxElements = {'Box', 'Median', 'Whisker', 'Upper Adjacent Value', 'Lower Adjacent Value', 'Upper Whisker', 'Lower Whisker'};
for i = 1:length(boxElements)
    boxHandles = findobj(h, 'Tag', boxElements{i});
    for j = 1:length(boxHandles)
        set(boxHandles(j), 'LineWidth', 2);
    end
end

% Make outliers more visible
outlierHandles = findobj(h, 'Tag', 'Outliers');
for j = 1:length(outlierHandles)
    set(outlierHandles(j), 'MarkerSize', 8, 'LineWidth', 1.5);
end

% Calculate appropriate y-axis limits that include the bounds
yMax = max(cellfun(@max, allData));
yMin = min(cellfun(@min, allData));
yRange = yMax - yMin;
extraSpace = yRange * 0.15;  % 15% extra space for text and markers
plotYMin = min(yMin - extraSpace, lower_bound - extraSpace/2);  % Ensure lower bound is visible
plotYMax = max(yMax + extraSpace/2, upper_bound + extraSpace/2);  % Ensure upper bound is visible

% Add grid and reference lines
%grid on
yl = yline([1/3], '-.r', 'β = 1/3', 'LineWidth', 2.5, 'FontSize', 14, 'FontWeight', 'bold');  % Changed from -1/3
yl.LabelHorizontalAlignment = 'left';
yl.LabelVerticalAlignment = 'bottom';

% Add equivalence bounds
yl_lower = yline([lower_bound], ':k', ['Lower bound: ' num2str(lower_bound, '%.3f')], 'LineWidth', 1.5, 'FontSize', 12);
yl_upper = yline([upper_bound], ':k', ['Upper bound: ' num2str(upper_bound, '%.3f')], 'LineWidth', 1.5, 'FontSize', 12);
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
             'FontWeight', 'bold', 'Color', [0 0.5 0], 'FontSize', 13);
    else
        text(i, starPos, 'non-equiv', 'HorizontalAlignment', 'center', ...
             'FontWeight', 'bold', 'Color', [0.8 0 0], 'FontSize', 13);
    end
end

% Restore grid and adjust layer order
grid off
set(gca, 'Layer', 'bottom')

% Set y-axis limits to include bounds and markers
ylim([plotYMin plotYMax])

% Titles and labels
title({'Comparison of Power Law β Values Across Analysis Methods', ...
       'Using Empirical Ellipse Tracing Data (Zarandi et al. 2023)'})

ylabel('Power Law β Value')
xlabel({' ', ...
    'Filter Types: BW (Butterworth), SG (Savitzky-Golay)', ...
    'Regression Types: LR (Linear), LMLS (Levenberg-Marquardt), IRLS (Iteratively Reweighted)', ...
    ' '})

% Add box plot explanation
boxExplain = {'Box Plot Elements:', ...
     '• Box: 25th to 75th percentiles', ...
     '• Center line: median', ...
     '• Notches: 95% CI of median', ...
     '• Whiskers: data range (excluding outliers)', ...
     '• Points: outliers'};

annotation('textbox', [.15 .7 .2 .2], 'String', boxExplain, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', ...
    'FontSize', 11, ...
    'FontWeight', 'bold', ... 
    'LineWidth', 1.5, ...
    'EdgeColor', [0.7 0.7 0.7]);

% Add equivalence testing explanation
eqvExplain = {'Equivalence Testing:', ...
     'equiv: mean within ±10% of 1/3', ...
     'non-equiv: mean outside ±10% bounds'};

annotation('textbox', [.7 .7 .2 .15], 'String', eqvExplain, ...
    'FitBoxToText', 'on', ...
    'BackgroundColor', 'white', ...
    'FontSize', 11, ...
    'FontWeight', 'bold', ... 
    'LineWidth', 1.5, ...
    'EdgeColor', [0.7 0.7 0.7]);

% Preserve title settings BEFORE calling setPublicationStyle
figTitle = get(fig, 'Name');  % Store the current title
titleHandle = title(gca, {'Comparison of Power Law β Values Across Analysis Methods', ...
            'Using Empirical Ellipse Tracing Data (Zarandi et al. 2023)'}, ...
            'FontSize', 16, ...
            'FontWeight', 'bold');

% Apply publication styling
setPublicationStyle(fig);

% Export with standardized settings
saveImage = true;  % Set to true to trigger export
if saveImage
    % Enable manual adjustment before export
    exportPublicationFigure(fig, figFilename11, ...
        'Width', 24, 'Height', 18, ...
        'Formats', {'png', 'eps'}, ...
        'Resolution', 1200, ...
        'ManualAdjust', true);  % Add manual adjustment parameter
end
