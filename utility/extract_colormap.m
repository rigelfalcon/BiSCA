function CM = extract_colormap(imagePath, numColor)
    % Extracts a colormap from an image with interpolation to a specified number of colors.
    %
    % Parameters:
    % imagePath: string, the path to the image file.
    % colorNum: integer, the number of colors to extract (optional, default is 256).
    %
    % Returns:
    % CM: matrix, the interpolated colormap.

    % Set a default colorNum if not provided

    if nargin < 1 || isempty(imagePath)
        error('extract_colormap:missingInput', ...
              'imagePath is required. Provide the path to a colorbar image.');
    end
    if nargin < 2 || isempty(numColor)
        numColor = 256; % Default to 256 colors, like turbo or hot in MATLAB
    end

    ext = get_ext(imagePath);
    % Check if the image is WebP format
    if strcmpi(ext, '.webp')
        oriPic = wpread(imagePath); % Use wpread for WebP files
    else
        oriPic = imread(imagePath); % Use imread for other image types
    end
    % Use rgb2ind to extract the colors
    [~, C] = rgb2ind(oriPic, 8);

    % Call the NTraveler function for color processing
    C = NTraveler(C);

    % Interpolate the colormap to 256 colors
    CM = interpColor(C, numColor);

end

function [ext] = get_ext(filename)
    [~, ~, ext] = fileparts(filename);
end