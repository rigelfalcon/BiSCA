function reflectiveScatter3(x, y, z, mksize, colors)
% reflectiveScatter3  Plot reflective points as small spheres.

% Check if colors are provided; if not, default to a single color
if nargin < 5 || isempty(colors)
    colors = repmat([0.5, 0.5, 0.5], length(x), 1); % Default gray color
end

% Create a sphere with 10 subdivisions for smoother visualization
[sx, sy, sz] = sphere(10);
mksize = mksize ./ 50000;

% Plot each point as a scaled and translated sphere
hold on;

for i = 1:length(x)
    surf(x(i) + sx * mksize, y(i) + sy * mksize, z(i) + sz * mksize, ...
        'FaceColor', colors(i, :), 'EdgeColor', 'none');
end

% Set lighting and material properties
light;
material metal; % Adjust to 'shiny' or 'dull' as needed
lighting gouraud; % Use Gouraud shading for smooth lighting effects
camlight headlight; % Add a camera light for better visualization

% Final touch to make the scene look better
axis equal;
hold off;
end
