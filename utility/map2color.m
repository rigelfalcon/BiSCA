function [colors,color_indices] = map2color(data, varargin)
%MAPTOCOLORMAP Maps a double vector to RGB colors based on a colormap.
%
%   COLORS = MAPTOCOLORMAP(DATA) maps the elements of the double vector
%   DATA to RGB colors using the default 'parula' colormap with 256 colors.
%
%   COLORS = MAPTOCOLORMAP(DATA, 'Colormap', CMAP_NAME) allows specifying
%   a different colormap. CMAP_NAME can be any valid MATLAB colormap name
%   (e.g., 'jet', 'hot', 'viridis', etc.).
%
%   COLORS = MAPTOCOLORMAP(DATA, 'NumColors', N) allows specifying the
%   number of colors in the colormap. Default is 256.
%
%   COLORS = MAPTOCOLORMAP(DATA, 'DataRange', [MIN, MAX]) allows specifying
%   the normalization range for DATA. Elements below MIN are mapped to the
%   first color, and elements above MAX are mapped to the last color.
%
%   Example:
%       data = rand(1, 100);
%       colors = mapToColormap(data, 'Colormap', 'jet', 'NumColors', 512);
%       scatter(1:100, rand(1, 100), 36, colors, 'filled');
%       colorbar;
%       colormap('jet');
%
%   See also COLORMAP, SCATTER, PLOT, IMAGESC.

    % Input parsing with name-value pairs
    p = inputParser;
    addRequired(p, 'data', @(x) isnumeric(x) && isvector(x));
    addParameter(p, 'Colormap', 'parula', @(x) ischar(x) || isstring(x));
    addParameter(p, 'NumColors', 256, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(p, 'DataRange', [], @(x) isempty(x) || (isnumeric(x) && numel(x) == 2));
    
    parse(p, data, varargin{:});
    
    data = p.Results.data(:); % Ensure data is a column vector
    cmapName = char(p.Results.Colormap);
    numColors = p.Results.NumColors;
    dataRange = p.Results.DataRange;
    
    % Select colormap
    try
        cmap = feval(cmapName, numColors);
    catch
        error('Invalid colormap name: %s', cmapName);
    end
    
    % Determine data range for normalization
    if isempty(dataRange)
        data_min = min(data);
        data_max = max(data);
    else
        data_min = dataRange(1);
        data_max = dataRange(2);
        if data_min >= data_max
            error('Invalid DataRange: Minimum value must be less than maximum value.');
        end
    end
    
    % Normalize data to [0, 1]
    normalized_data = (data - data_min) / (data_max - data_min);
    
    % Clip normalized data to [0, 1]
    normalized_data = max(min(normalized_data, 1), 0);
    
    % Map normalized data to colormap indices
    color_indices = round(normalized_data * (numColors - 1)) + 1;
    
    % Assign RGB colors
    colors = cmap(color_indices, :);
end
