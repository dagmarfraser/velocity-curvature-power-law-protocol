function [colors, lineStyles, markerStyles] = getConsistentStyles()
    % Return consistent colors, line styles, and marker styles for all figures
    %
    % Created for the velocity-curvature power law protocol
    % This function provides standardized visual elements across all figures
    % 
    % Returns:
    % --------
    % colors : matrix
    %     Matrix of RGB colors (n x 3 matrix)
    % lineStyles : cell array
    %     Cell array of line style strings
    % markerStyles : cell array
    %     Cell array of marker style strings
    
    % Define a colorblind-friendly color palette
    colors = [
        0.0000, 0.4470, 0.7410;  % Blue
        0.8500, 0.3250, 0.0980;  % Orange
        0.4940, 0.1840, 0.5560;  % Purple
        0.4660, 0.6740, 0.1880;  % Green
        0.3010, 0.7450, 0.9330;  % Cyan
        0.6350, 0.0780, 0.1840;  % Dark red
        0.9290, 0.6940, 0.1250   % Yellow
    ];
    
    % Define line styles
    lineStyles = {'-', '--', '-.', ':'};
    
    % Define marker styles
    markerStyles = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
end