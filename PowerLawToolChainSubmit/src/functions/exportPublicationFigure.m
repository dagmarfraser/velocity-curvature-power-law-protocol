function exportPublicationFigure(fig, filename, varargin)
    % Export figures in publication-quality formats
    %
    % Created for the velocity-curvature power law protocol
    % This function standardizes figure export for:
    % - Manuscript submission 
    % - Journal publication in Exp Brain Research
    %
    % Parameters:
    % -----------
    % fig : figure handle
    %     The figure to export
    % filename : string
    %     Base filename (without extension)
    % varargin : parameter/value pairs
    %     Optional parameters:
    %     - 'Resolution': Output resolution in DPI (default: 1200)
    %     - 'Width': Figure width in centimeters (default: 16)
    %     - 'Height': Figure height in centimeters (default: 12)
    %     - 'Formats': Cell array of formats to export (default: {'png', 'eps'})
    %     - 'Directory': Directory to save figures (default: '[project_root]/figures')
    %     - 'ManualAdjust': Enable manual adjustment mode before export (default: false)
    
    % Parse optional input arguments
    p = inputParser;
    addParameter(p, 'Resolution', 1200, @isnumeric);
    addParameter(p, 'Width', 16, @isnumeric);
    addParameter(p, 'Height', 12, @isnumeric);
    addParameter(p, 'Formats', {'png', 'eps'}, @iscell);
    addParameter(p, 'ManualAdjust', false, @islogical);
    
    % Get default figures directory (parent's parent of current directory/figures)
    currentDir = pwd;
    % Navigate up to PowerLawToolChainSubmit directory
    if contains(currentDir, 'src')
        % If we're in the src directory
        parentDir = fileparts(currentDir);
    else
        % If we're already in PowerLawToolChainSubmit or somewhere else
        parentDir = currentDir;
    end
    defaultFiguresDir = fullfile(parentDir, 'figures');
    
    addParameter(p, 'Directory', defaultFiguresDir, @ischar);
    parse(p, varargin{:});
    
    % Get parameters
    res = p.Results.Resolution;
    width = p.Results.Width;
    height = p.Results.Height;
    formats = p.Results.Formats;
    folderSave = p.Results.Directory;
    manualAdjust = p.Results.ManualAdjust;
    
    % Ensure figure has correct styling
    %setPublicationStyle(fig);
    
    % Set figure size
    set(fig, 'Units', 'centimeters', 'Position', [1, 1, width, height], 'PaperPosition', [0, 0, width, height]);
    
    % Get parameters
    res = p.Results.Resolution;
    width = p.Results.Width;
    height = p.Results.Height;
    formats = p.Results.Formats;
    folderSave = p.Results.Directory;
    manualAdjust = p.Results.ManualAdjust;
    
    % Create directory if it doesn't exist
    if ~exist(folderSave, 'dir')
        mkdir(folderSave);
    end
    
    % Enable manual adjustment if requested
    if manualAdjust
        plotedit(fig, 'on');
        disp('PRESS ANY KEY ONCE PLOT ELEMENTS MOVED TO FINAL POSITION!');
        pause
        plotedit(fig, 'off');
    end
    
    % Prepare figure for export
    set(fig, 'PaperPositionMode', 'auto'); 
    set(fig, 'InvertHardcopy', 'off');  % Preserve background color
    
    % Make sure tick rotation is preserved
    axesHandles = findall(fig, 'Type', 'axes');
    for i = 1:length(axesHandles)
        ax = axesHandles(i);
        % Skip invisible or special axes (like colorbars)
        if strcmpi(get(ax, 'Visible'), 'off') || ...
                strcmp(get(ax, 'Tag'), 'legend') || ...
                strcmp(get(ax, 'Tag'), 'colorbar')
            continue;
        end
        % Store rotations before export
        xrot = get(ax, 'XTickLabelRotation');
        yrot = get(ax, 'YTickLabelRotation');
        setappdata(ax, 'XTickLabelRotation', xrot);
        setappdata(ax, 'YTickLabelRotation', yrot);
    end
    
    % Export in specified formats
    for i = 1:length(formats)
        format = lower(formats{i});
        outputFile = fullfile(folderSave, [filename, '.', format]);
        
        % Use appropriate export method based on format
        try
            % Restore the tick rotations immediately before export
            for j = 1:length(axesHandles)
                ax = axesHandles(j);
                % Skip special axes
                if strcmpi(get(ax, 'Visible'), 'off') || ...
                        strcmp(get(ax, 'Tag'), 'legend') || ...
                        strcmp(get(ax, 'Tag'), 'colorbar')
                    continue;
                end
                % Restore rotations
                if isappdata(ax, 'XTickLabelRotation')
                    xrot = getappdata(ax, 'XTickLabelRotation');
                    set(ax, 'XTickLabelRotation', xrot);
                end
                if isappdata(ax, 'YTickLabelRotation')
                    yrot = getappdata(ax, 'YTickLabelRotation');
                    set(ax, 'YTickLabelRotation', yrot);
                end
            end
            
            switch format
                case 'png'
                    % Using ContentType vector preserves vector elements in the PNG where possible
                    exportgraphics(fig, outputFile, ...
                        'Resolution', res, ...
                        'BackgroundColor', 'white', ...
                        'ContentType', 'vector');
                    
                case 'eps'
                    exportgraphics(fig, outputFile, ...
                        'Resolution', res, ...
                        'BackgroundColor', 'white');
                    
                case 'pdf'
                    exportgraphics(fig, outputFile, ...
                        'Resolution', res, ...
                        'BackgroundColor', 'white', ...
                        'ContentType', 'vector');
                    
                case 'svg'
                    exportgraphics(fig, outputFile, ...
                        'ContentType', 'vector');
                    
                case 'tiff'
                    exportgraphics(fig, outputFile, ...
                        'Resolution', res, ...
                        'BackgroundColor', 'white');
            end
            fprintf('Saved: %s\n', outputFile);
        catch ME
            warning('Error exporting %s: %s', outputFile, ME.message);
            
            % Fall back to print function for older MATLAB versions
            if contains(ME.message, 'exportgraphics')
                try
                    if strcmp(format, 'png')
                        print(fig, outputFile, '-dpng', sprintf('-r%d', res));
                    elseif strcmp(format, 'eps')
                        print(fig, outputFile, '-depsc', sprintf('-r%d', res));
                    elseif strcmp(format, 'pdf')
                        print(fig, outputFile, '-dpdf', sprintf('-r%d', res));
                    elseif strcmp(format, 'tiff')
                        print(fig, outputFile, '-dtiff', sprintf('-r%d', res));
                    end
                    fprintf('Saved using print function: %s\n', outputFile);
                catch ME2
                    warning('Fallback export also failed: %s', ME2.message);
                end
            end
        end
    end
    
    % Display output location
    fprintf('Figure(s) saved to: %s\n', folderSave);
end