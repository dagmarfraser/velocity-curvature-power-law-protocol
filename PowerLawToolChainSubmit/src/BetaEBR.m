% Create a figure with larger dimensions and better visibility
figure('Position', [100 100 1000 600]);

% Angular frequencies from ShapesHuhEBR.m
nuList = [0.01 2/5 3/5 2/3 4/5 4/3 3/2 2 5/2 3 4 6];
nuListCell = {'0.01' '2/5' '3/5' '2/3' '4/5' '4/3' '3/2' '2' '5/2' '3' '4' '6'};

% Create a denser list for the continuous line, using logarithmic spacing
nu_dense = logspace(-2, log10(6), 1000);

% Calculate beta for each nu using the formula
beta = @(nu) (2/3) * ((1 + nu.^2/2) ./ (1 + nu.^2 + nu.^4/15));

beta_dense = beta(nu_dense);
beta_points = beta(nuList);

% Create the plot with both line and points
hold on
% Plot continuous line
plot(nu_dense, beta_dense, 'LineWidth', 2, 'Color', [0.2 0.6 0.8]);
% Plot specific points
scatter(nuList, beta_points, 100, 'ro', 'filled', 'MarkerEdgeColor', 'k');

% Add horizontal line at β = 1/3
yline(1/3, '--k', '\beta = 1/3', 'LineWidth', 1.5, 'FontSize', 14, 'LabelHorizontalAlignment', 'left');

% Customize the plot
grid on
box on
set(gca, 'XScale', 'log')  % Set x-axis to logarithmic scale
xlabel('Angular Frequency (\nu)', 'FontSize', 14, 'Interpreter', 'tex');
ylabel('Power Law Exponent (\beta)', 'FontSize', 14, 'Interpreter', 'tex');
title('Power Law Exponent vs Angular Frequency', 'FontSize', 16);

% Add annotations for specific points with adjusted positions for log scale
for i = 1:length(nuList)
    % Add some vertical offset to prevent overlap
    text(nuList(i), beta_points(i) + 0.02, nuListCell{i}, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 12);
end

% Set axis limits with some padding
axLimits = axis;
set(gca, 'XLim', [0.008 7], 'YLim', [0.1 0.9]); % Further expanded y-axis limits

% Add legend
legend('Continuous \beta(\nu)', 'Notable \nu values', 'Location', 'northeast', 'FontSize', 12);

% Add explanatory text with improved LaTeX formatting - much larger and bold
text(0.015, 0.8, {'$\mathbf{\beta(\nu) = \frac{2}{3} \cdot \frac{1 + \frac{\nu^2}{2}}{1 + \nu^2 + \frac{\nu^4}{15}}}$'}, ...
    'FontSize', 28, 'Interpreter', 'latex', 'HorizontalAlignment', 'left');

% Customize grid for log scale
set(gca, 'XGrid', 'off', 'YGrid', 'on')  % Turn off default x grid
grid minor  % Add minor grid lines

% Add custom major grid lines at specific ν values
axisLims = axis;
for nu = [0.01 0.1 1 10]
    if nu >= axisLims(1) && nu <= axisLims(2)
        line([nu nu], [axisLims(3) axisLims(4)], 'Color', [0.8 0.8 0.8], 'LineStyle', '-', 'LineWidth', 0.5);
    end
end

hold off