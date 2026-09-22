function figHandle = boxscatter(data, markerStyle, markerColor, markerSize, lineWidth)
%BOXSCATTER Boxplot with overlaid jittered scatter points.
%
% DESCRIPTION:
%   BOXPLOT + SCATTER visualization for grouped data. Each column of the
%   input matrix is treated as a group. A standard boxplot is drawn, and
%   individual data points are overlaid using jittered scatter markers.
%   Marker appearance (style, color, size, transparency) and boxplot line
%   properties can be customized.
%
% SYNTAX:
%   figHandle = boxscatter(data)
%   figHandle = boxscatter(data, markerStyle)
%   figHandle = boxscatter(data, markerStyle, markerColor)
%   figHandle = boxscatter(data, markerStyle, markerColor, markerSize)
%   figHandle = boxscatter(data, markerStyle, markerColor, markerSize, lineWidth)
%
% INPUTS:
%   data         - [nSamples x nGroups] numeric matrix.
%                  Each column corresponds to one group.
%
%   markerStyle  - (optional) Marker symbol(s) for scatter points.
%                  Example: 'o', 's', '^'. Can be a vector of styles.
%                  Default: 'o'
%
%   markerColor  - (optional) [nGroups x 3] RGB matrix specifying colors.
%                  Default: lines(nGroups)
%
%   markerSize   - (optional) Scalar specifying marker size.
%                  Default: 36
%
%   lineWidth    - (optional) Scalar specifying line width of boxplot elements.
%                  Default: 0.5
%
% OUTPUT:
%   figHandle    - Handle to the generated figure.
%
% FEATURES:
%   - Jittered scatter to reduce overlap
%   - Custom marker styles per group
%   - Matching outlier colors to box colors
%   - Adjustable marker transparency (hardcoded in function)
%
% EXAMPLE:
%   % Generate sample data (3 groups)
%   data = randn(30,3) + [0 1 2];
%
%   % Define colors
%   colors = [0.2 0.4 0.8;
%             0.8 0.3 0.3;
%             0.2 0.7 0.4];
%
%   % Plot
%   boxscatter(data, 'o', colors, 40, 1);
%
%   % Result: Boxplot with colored, semi-transparent scatter points
%
% AUTHOR: Edgar Guevara, PhD
% Custom function for enhanced boxplot visualization.

    % Check for optional arguments and set defaults
    if nargin < 2 || isempty(markerStyle)
        markerStyle = 'o'; % Default marker style
    end
    if nargin < 3 || isempty(markerColor)
        markerColor = lines(size(data, 2)); % Default colors
    end
    if nargin < 4 || isempty(markerSize)
        markerSize = 36; % Default marker size
    end
    if nargin < 5 || isempty(lineWidth)
        lineWidth = 0.5; % Default line width for boxplot
    end
    
    % Transparency value for scatter markers
    % 1 → fully opaque; 0 → fully transparent
    alphaVal = 0.5;

    % Create the boxplot
    boxplotHandle = boxplot(data);
    figHandle = gcf;
    ax=gca;

    % Change median lines to dotted lines
    % Find all median line objects
    medianLines = findobj(boxplotHandle, 'Tag', 'Median');
    % Change the line style of median lines to dotted
    set(medianLines, 'LineStyle', '--', 'Color', 'k', 'LineWidth', lineWidth); % 'k' for black color, '--' for dashed line style

     % Set line width for boxplot
    h = findobj(boxplotHandle, 'Tag', 'Box');
    for idx = 1:length(h)
        set(h(idx), 'LineWidth', lineWidth, 'Color', 'k');
    end
    % Set line color for median line
    % hMedian = findobj(gca, 'Tag', 'Median');
    % for idx = length(hMedian):-1:1
    %     set(hMedian(idx), 'Color', markerColor(idx,:),'LineWidth', lineWidth);
    % end
    % Set line width for whiskers
    hWhisker = findobj(ax, 'Tag', 'Whisker');
    set(hWhisker, 'LineWidth', lineWidth);
    
    % Set outlier colors to match each box
    hOutliers = findobj(ax, 'Tag', 'Outliers');

    % Outliers are returned in reverse order → flip for correct mapping
    hOutliers = flipud(hOutliers);

    for idx = 1:min(length(hOutliers), size(markerColor,1))
        set(hOutliers(idx), ...
            'MarkerEdgeColor', markerColor(idx,:), ...
            'MarkerFaceColor', markerColor(idx,:)); % optional (filled markers)
    end

    hold on; % Keep the boxplot visible when adding the scatter plot
    
    % Generate random x-coordinates for scatter plot
    groups = size(data, 2);
    xCoords = repmat(1:groups, size(data, 1), 1) + (rand(size(data)) - 0.5) * 0.4;
    
    % Plot the scatter plot for each group
    for idx = 1:groups
        % Determine the current marker style
        if length(markerStyle) >= idx
            currentMarkerStyle = markerStyle(idx);
        else
            currentMarkerStyle = markerStyle(end); % Use the last one if not enough styles provided
        end
        
        % Plot scatter for current group
        scatter(xCoords(:, idx), data(:, idx), 'SizeData', markerSize, ...
            'Cdata', markerColor(idx, :), 'Marker', currentMarkerStyle, ...
            'LineWidth', lineWidth, 'MarkerFaceAlpha', alphaVal, 'MarkerEdgeAlpha', alphaVal);
        % 1 → fully opaque; 0 → fully transparent
    end

   
    
    hold off; % Release the plot for further commands
end

% EOF