function setPublicationStyle(fig)
% Apply publication-quality styling to all elements of a figure
%
% Created for the velocity-curvature power law protocol
% This function standardizes the appearance of all figures for:
% - Manuscript submission
% - Journal publication in Exp Brain Research
%
% Parameters:
% -----------
% fig : figure handle
%     The figure to apply styling to

% Base styling parameters
fontName = 'Arial';
baseFontSize = 14;
titleFontSize = 16;
lineWidth = 2.0;
markerSize = 8;
axisLineWidth = 1.5;

% Set figure properties
set(fig, 'Color', 'white', ...
    'PaperPositionMode', 'auto');

% Find all axes in the figure (including subplots)
axesHandles = findall(fig, 'Type', 'axes');

% Capture existing titles before styling
titleHandles = [];
titleStrings = {};
titleProps = {};

for i = 1:length(axesHandles)
    ax = axesHandles(i);
    titleHandle = get(ax, 'Title');

    if ~isempty(titleHandle)
        % Store title handles and properties
        titleHandles(end+1) = titleHandle;
        titleStrings{end+1} = get(titleHandle, 'String');

        % Capture key properties to restore
        titleProps{end+1} = struct(...
            'String', get(titleHandle, 'String'), ...
            'FontSize', get(titleHandle, 'FontSize'), ...
            'FontWeight', get(titleHandle, 'FontWeight'), ...
            'Color', get(titleHandle, 'Color'), ...
            'Interpreter', get(titleHandle, 'Interpreter'));
    end
end

% Run the rest of the publication styling function
for i = 1:length(axesHandles)
    ax = axesHandles(i);

    % Skip invisible or special axes (like colorbars)
    if strcmpi(get(ax, 'Visible'), 'off') || ...
            strcmp(get(ax, 'Tag'), 'legend') || ...
            strcmp(get(ax, 'Tag'), 'colorbar')
        continue;
    end

    % Set axis properties
    set(ax, 'FontName', fontName, ...
        'FontSize', baseFontSize, ...
        'LineWidth', axisLineWidth, ...
        'Box', 'on', ...
        'TickDir', 'out', ...
        'TitleFontSizeMultiplier', 1.2, ...
        'TickLength', [0.015 0.025], ...
        'TitleFontWeight', 'bold', ...
        'XGrid', 'off', ...
        'YGrid', 'off', ...
        'XAxisLocation', 'bottom', ...
        'YAxisLocation', 'left');
end

% Restore titles after styling
for i = 1:length(titleHandles)
    try
        % Recreate title with preserved properties
        currTitle = titleProps{i};

        % Use title() to ensure the title is properly attached
        titleH = title(axesHandles(i), currTitle.String, ...
            'FontName', fontName, ...
            'FontSize', titleFontSize, ...  % Use predefined title font size
            'FontWeight', currTitle.FontWeight, ...
            'Color', currTitle.Color, ...
            'Interpreter', currTitle.Interpreter);
    catch ME
        warning('Failed to restore title: %s', ME.message);
    end
end

% Rest of the existing styling function remains the same
% (line styling, scatter plots, patches, etc.)

% [Keep the rest of the existing function content here]
% ... (legend styling, suptitle styling, etc.)


% Ensure no top axis ticks
if isprop(ax, 'XRuler')
    try
        ax.XRuler.TicksMode = 'manual';
        ax.XRuler.TicksRequestedPoints = get(ax, 'XTick');
        ax.XRuler.Axle.Visible = 'off';
    catch
        % Fallback for older MATLAB versions
        set(ax, 'XTickMode', 'manual');
    end
end

% Remove top tick marks to avoid intersecting with titles
% This approach works in most MATLAB versions

% Get current tick positions
xticks = get(ax, 'XTick');
yticks = get(ax, 'YTick');

% Store original axes position
originalUnits = get(ax, 'Units');
set(ax, 'Units', 'normalized');
axPos = get(ax, 'Position');

% Increase spacing between plot and title
% This adds extra space to prevent ticks from overlapping with title
% Adjust title spacing more carefully
if ~isempty(get(ax, 'Title'))
    titleH = get(ax, 'Title');
    titleStr = get(titleH, 'String');
    if ~isempty(titleStr)
        % Get current axes position
        originalUnits = get(ax, 'Units');
        set(ax, 'Units', 'normalized');
        axPos = get(ax, 'Position');

        % Reduce height slightly but less aggressively
        axPos(4) = axPos(4) * 0.97;  % Reduce by only 3% instead of 5%
        set(ax, 'Position', axPos);

        % % Adjust title position more subtly
        % titlePos = get(titleH, 'Position');
        % if length(titlePos) >= 2
        %     titlePos(2) = titlePos(2) * 1.02;  % Move up by only 2%
        %     set(titleH, 'Position', titlePos);
        % end

        % Restore original units
        set(ax, 'Units', originalUnits);
    end
end

% Restore original units
set(ax, 'Units', originalUnits);

% Use a simpler approach for removing top ticks
% Just set the tick direction and adjust position slightly
set(ax, 'TickDir', 'in');

% Simple way to create more space between title and plot
% Adjust the top margin of the plot slightly
ylim = get(ax, 'YLim');
ydiff = diff(ylim);
set(ax, 'YLim', [ylim(1), ylim(2) + ydiff * 0.05])

% Preserve x-tick label rotation if it's already set
currentXRotation = get(ax, 'XTickLabelRotation');
if ~isempty(currentXRotation) && currentXRotation ~= 0
    % Remember the rotation value to reapply after other styling
    xRotationValue = currentXRotation;
else
    xRotationValue = 0; % Default is horizontal
end

% Set title properties if it exists
titleHandle = get(ax, 'Title');
if ~isempty(titleHandle)
    set(titleHandle, 'FontSize', titleFontSize, ...
        'FontName', fontName, ...
        'FontWeight', 'bold');

    % Add space between title and plot
    % Get and increase vertical position for title spacing
    titlePos = get(titleHandle, 'Position');
    if ~isempty(titlePos) && length(titlePos) >= 3
        % Only adjust non-empty titles
        titleStr = get(titleHandle, 'String');
        % if ~isempty(titleStr)
        %     % Increase vertical position by 5%
        %     titlePos(2) = titlePos(2) * 1.05;
        %     set(titleHandle, 'Position', titlePos);
        % end
    end
end

% Set label properties
set(get(ax, 'XLabel'), 'FontSize', baseFontSize, 'FontName', fontName);
set(get(ax, 'YLabel'), 'FontSize', baseFontSize, 'FontName', fontName);

% Set all lines to have consistent width
lines = findall(ax, 'Type', 'line');
for j = 1:length(lines)
    % Skip reference lines like axis borders
    if ~strcmp(get(lines(j), 'Tag'), '')
        continue;
    end

    set(lines(j), 'LineWidth', lineWidth);

    % Set marker size if marker exists
    if ~isempty(get(lines(j), 'Marker')) && ~strcmpi(get(lines(j), 'Marker'), 'none')
        set(lines(j), 'MarkerSize', markerSize);
    end
end

% Style scatter plots
scatters = findall(ax, 'Type', 'scatter');
for j = 1:length(scatters)
    set(scatters(j), 'LineWidth', lineWidth);
    set(scatters(j), 'SizeData', markerSize^2 * 4);
end

% Set all patch elements (for filled areas, boxplots, etc.)
patches = findall(ax, 'Type', 'patch');
for j = 1:length(patches)
    set(patches(j), 'LineWidth', lineWidth/2);
end

% Handle boxplot elements separately for better publication quality
boxElements = {'Box', 'Median', 'Whisker', 'Upper Adjacent Value', 'Lower Adjacent Value', 'Upper Whisker', 'Lower Whisker'};
for j = 1:length(boxElements)
    boxHandles = findobj(ax, 'Tag', boxElements{j});
    if ~isempty(boxHandles)
        for k = 1:length(boxHandles)
            set(boxHandles(k), 'LineWidth', lineWidth);
        end
    end
end

% Make boxplot outliers more visible
outlierHandles = findobj(ax, 'Tag', 'Outliers');
if ~isempty(outlierHandles)
    for k = 1:length(outlierHandles)
        set(outlierHandles(k), 'MarkerSize', markerSize, 'LineWidth', lineWidth/1.5);
    end
end

% Set all text elements to use the same font
textElements = findall(ax, 'Type', 'text');
for j = 1:length(textElements)
    % Don't change font size for specific annotations that might need custom sizing
    currentFontSize = get(textElements(j), 'FontSize');
    set(textElements(j), 'FontName', fontName);
end

% Add padding to axis limits to prevent content crowding
% Only apply to non-categorical axes
if ~isa(ax.XAxis, 'matlab.graphics.axis.CategoricalRuler')
    xlims = get(ax, 'XLim');
    xrange = diff(xlims);
    set(ax, 'XLim', [xlims(1) - xrange*0.02, xlims(2) + xrange*0.02]);
end

if ~isa(ax.YAxis, 'matlab.graphics.axis.CategoricalRuler')
    ylims = get(ax, 'YLim');
    yrange = diff(ylims);
    set(ax, 'YLim', [ylims(1) - yrange*0.05, ylims(2) + yrange*0.05]);
end

% Reapply x-tick rotation if it was set
if xRotationValue ~= 0
    set(ax, 'XTickLabelRotation', xRotationValue);
end


% Set legend properties if it exists
legendHandles = findobj(fig, 'Type', 'legend');
for i = 1:length(legendHandles)
    set(legendHandles(i), 'FontSize', baseFontSize, ...
        'FontName', fontName, ...
        'Box', 'on', ...
        'EdgeColor', [0.7 0.7 0.7]);

    % Make legend line widths consistent
    legLines = findobj(legendHandles(i), 'Type', 'line');
    for j = 1:length(legLines)
        set(legLines(j), 'LineWidth', lineWidth);
    end
end

% Set figure title (suptitle) if it exists
titleTexts = findall(fig, 'Type', 'text');
for i = 1:length(titleTexts)
    if strcmp(get(titleTexts(i), 'Tag'), 'suptitle')
        set(titleTexts(i), 'FontName', fontName, ...
            'FontSize', titleFontSize+2, ...
            'FontWeight', 'bold');
    end
end
