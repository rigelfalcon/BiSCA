function varargout = surfscatter(faces, vertices, vin, node, opt)
% SURFACEPLOT  Plot brain surface with optional node data for both hemispheres
%   Inputs:
%       faces     - Cell array of face matrices for left and right hemispheres
%       vertices  - Cell array of vertex coordinates for left and right hemispheres
%       vin       - Vertex or node values (single array or cell array)
%       node      - Optional node coordinates
%       opt       - Options structure with fields:
%                   - hemisphere: 'both' (default), 'left', or 'right'
%                   - reflective: Boolean for reflective scatter (default: false)
%                   - isbinary: Boolean for binary colormap (default: false)
%                   - cmp: Custom colormap (optional)

% Set default options
opt = set_defaults(opt, 'hemisphere', 'both', 'reflective', false, 'isbinary', false);
opt = set_defaults(opt, 'std_range', false);
opt = set_defaults(opt, 'tess', []);

% Load default brain data if vertices or faces are empty
if isempty(vertices) || isempty(faces)
    % Prefer opt.tess (struct with Vertices, Faces, Atlas).
    if ~isempty(opt.tess)
        tess = opt.tess;
    else
        repo_root = fileparts(fileparts(mfilename('fullpath')));
        mesh_path = fullfile(repo_root, 'data', 'sEEG', 'anatomy', 'anat', ...
            'tess_cortex_mid_low.mat');
        if exist(mesh_path, 'file') ~= 2
            error('surfscatter:noMesh', ...
                'No cortex mesh available. Pass opt.tess or place tess_cortex_mid_low.mat at:\n  %s', ...
                mesh_path);
        end
        tess = load(mesh_path);
    end

    vertices{1} = tess.Vertices(tess.Atlas(2).Scouts(1).Vertices, :);
    vertices{2} = tess.Vertices(tess.Atlas(2).Scouts(2).Vertices, :);

    left_vertex_indices_new = 1:length(tess.Atlas(2).Scouts(1).Vertices);
    right_vertex_indices_new = 1:length(tess.Atlas(2).Scouts(2).Vertices);

    [~, loc] = ismember(tess.Faces, tess.Atlas(2).Scouts(1).Vertices);
    valid_entries = all(loc > 0, 2);
    faces{1} = left_vertex_indices_new(loc(valid_entries, :));

    [~, loc] = ismember(tess.Faces, tess.Atlas(2).Scouts(2).Vertices);
    valid_entries = all(loc > 0, 2);
    faces{2} = right_vertex_indices_new(loc(valid_entries, :));

    if opt.inflated
        opt.inflated = false; % Inflated mesh not bundled; fall back to non-inflated.
    end

end

% Handle case where vertices/faces are not cell arrays
if ~iscell(vertices)
    [faces, vertices] = separate_hemispheres(faces, vertices);
end

% Compute vertex values from nodes if provided, otherwise use vin directly
if ~isempty(node)
    radius = opt.radius; % Assume radius is provided or defaulted elsewhere

    if iscell(vin)
        v{1} = node2vertice(vin{1}, vertices{1}, node{1}, radius);
        v{2} = node2vertice(vin{2}, vertices{2}, node{2}, radius);
    else
        v{1} = node2vertice(vin, vertices{1}, node, radius);
        v{2} = node2vertice(vin, vertices{2}, node, radius);
    end

else

    if iscell(vin)
        v = vin;
    else
        v{1} = vin; % Assume vin applies to both if not separated
        v{2} = vin;
    end

end

% Plotting
% clf
% set(gcf, 'Position', [100, 100, 800, 600]);
mksize = 10;

if isfield(opt, 'mksize') && ~isempty(opt.mksize)
    mksize = opt.mksize;
end

ignoreoutlier = true; % Ignore outliers in surfpatch

hemisphere = opt.hemisphere;

if strcmp(hemisphere, 'both')
    % Plot both hemispheres
    p1 = surfpatch(faces{1}, vertices{1}, v{1}, ignoreoutlier);
    hold on
    p2 = surfpatch(faces{2}, vertices{2}, v{2}, ignoreoutlier);

    if ~isempty(node)

        if opt.reflective
            cin = get_colors(vin, cmp);
            reflectiveScatter3(node(:, 1), node(:, 2), node(:, 3), mksize, cin);
        else
            scatter3(node(:, 1), node(:, 2), node(:, 3), mksize, vin, 'filled', 'o');
        end

    end

    view(-90, 90); % Top view showing both hemispheres
elseif strcmp(hemisphere, 'left')
    % Plot left hemisphere only
    p = surfpatch(faces{1}, vertices{1}, v{1}, ignoreoutlier);

    if ~isempty(node)
        idx = node(:, 2) > 0; % Assuming y>0 is left hemisphere

        if opt.reflective
            cin = get_colors(vin(idx), cmp);
            reflectiveScatter3(node(idx, 1), node(idx, 2), node(idx, 3), mksize, cin);
        else
            scatter3(node(idx, 1), node(idx, 2), node(idx, 3), mksize, vin(idx), 'filled', 'o');
        end

    end

    view(-90, 0); % Lateral view of left hemisphere
elseif strcmp(hemisphere, 'right')
    % Plot right hemisphere only
    p = surfpatch(faces{2}, vertices{2}, v{2}, ignoreoutlier);

    if ~isempty(node)
        idx = node(:, 2) <= 0; % Assuming y<=0 is right hemisphere

        if opt.reflective
            cin = get_colors(vin(idx), cmp);
            reflectiveScatter3(node(idx, 1), node(idx, 2), node(idx, 3), mksize, cin);
        else
            scatter3(node(idx, 1), node(idx, 2), node(idx, 3), mksize, vin(idx), 'filled', 'o');
        end

    end

    view(90, 0); % Lateral view of right hemisphere
else
    error('Invalid hemisphere option. Choose ''both'', ''left'', or ''right''.');
end

% Set axis limits
if strcmp(hemisphere, 'both')

    if ~isempty(node)
        all_x = [vertices{1}(:, 1); vertices{2}(:, 1); node(:, 1)];
        all_y = [vertices{1}(:, 2); vertices{2}(:, 2); node(:, 2)];
        all_z = [vertices{1}(:, 3); vertices{2}(:, 3); node(:, 3)];
    else
        all_x = [vertices{1}(:, 1); vertices{2}(:, 1)];
        all_y = [vertices{1}(:, 2); vertices{2}(:, 2)];
        all_z = [vertices{1}(:, 3); vertices{2}(:, 3)];
    end

elseif strcmp(hemisphere, 'left')

    if ~isempty(node)
        idx = node(:, 2) > 0;
        all_x = [vertices{1}(:, 1); node(idx, 1)];
        all_y = [vertices{1}(:, 2); node(idx, 2)];
        all_z = [vertices{1}(:, 3); node(idx, 3)];
    else
        all_x = vertices{1}(:, 1);
        all_y = vertices{1}(:, 2);
        all_z = vertices{1}(:, 3);
    end

elseif strcmp(hemisphere, 'right')

    if ~isempty(node)
        idx = node(:, 2) <= 0;
        all_x = [vertices{2}(:, 1); node(idx, 1)];
        all_y = [vertices{2}(:, 2); node(idx, 2)];
        all_z = [vertices{2}(:, 3); node(idx, 3)];
    else
        all_x = vertices{2}(:, 1);
        all_y = vertices{2}(:, 2);
        all_z = vertices{2}(:, 3);
    end

end

xlimits = 1.2 * [min(all_x), max(all_x)];
ylimits = 1.2 * [min(all_y), max(all_y)];
zlimits = 1.2 * [min(all_z), max(all_z)];
xlim(xlimits);
ylim(ylimits);
zlim(zlimits);

% Set color limits and colormap
sigma = std(vin);
mu = mean(vin);
% climits = double([min(vin), max(vin)]);
climits = double([min([v{1};v{2}]), max([v{1};v{2}])]);

if sigma > 0 && ~opt.isbinary && opt.std_range
    climits = gather([mu - 1 * sigma, mu + 1 * sigma]);
end

% Guard against degenerate clim (e.g., all-zero data)
if climits(1) >= climits(2)
    climits(2) = climits(1) + 1;
end

cmp = plt.inferno;%redbluecmap;

if opt.isbinary
    cmp = flip(autumn(100), 1);
    cmp = [cmp(1, :); cmp(end, :)];
end

if isfield(opt, 'cmp') && ~isempty(opt.cmp)
    cmp = opt.cmp;
end

colormap(cmp);
clim(climits);

% Add colorbar
% c = colorbar;
% c.Location = 'southoutside';
% c.FontSize = 14;
% c = colorbar;
% c.Location = 'southoutside';
% c.FontSize = 14;
% pos = c.Position;
% original_center = pos(1) + pos(3) / 2;
% new_width = pos(3) * 0.5; % Make it half as wide
% new_left = original_center - new_width / 2;
% c.Position = [new_left, pos(2), new_width, pos(4)];

% Handle output arguments
if nargout == 1

    if strcmp(hemisphere, 'both')
        varargout{1} = [p1, p2];
    else
        varargout{1} = p;
    end

end

if nargout >= 2
    varargout{2} = v;
end

end

% Helper function to get colors
function cin = get_colors(vin, cmp)

if ischar(cmp) || isstring(cmp)
    cin = eval([cmp, '(', num2str(length(unique(vin))), ')']);
elseif isnumeric(cmp) && size(cmp, 2) == 3
    cin = cmp;
end

[~, cidx] = map2color(vin, 'NumColors', size(cin, 1));
cin = cin(cidx, :);
end

% Note: Assumes surfpatch, node2vertice, separate_hemispheres, reflectiveScatter3,
% redbluecmap, and map2color are defined elsewhere.
% Helper function to set default options
function opt = set_defaults(opt, varargin)

for i = 1:2:length(varargin)

    if ~isfield(opt, varargin{i}) || isempty(opt.(varargin{i}))
        opt.(varargin{i}) = varargin{i + 1};
    end

end

end

function v = node2vertice(vin, vertices, node, radius)
% v=NwSmooth(node,vin,radius,vertices);
if nargin < 4 || isempty(radius)
    radius = max(max(vertices, [], 1) - min(vertices, [], 1) / 2);
    % radius=radius*ones(size(vertices,2),1);
    m = 100;
    hlist = get_hList(m, [radius / 100, radius * 5], @logspace);%logspace
    % [~, nw_gcv] = NwSmoothGCV(gather(node), gather(vin), gather(hlist));
    [~, nw_gcv] = NwSmoothGCVOriginalScaleLogGaussian(gather(node), gather(vin), gather(hlist));
    % [nw_gcv]=cvNwSmooth(node,vin,node,m,hlist);
    v = NwSmooth(node, vin, nw_gcv.hmin, vertices);
else
    v = NwSmooth(node, vin, radius, vertices);
end

end

function p = surfpatch(faces, vertices, v, ignoreoutlier)
p = patch('Faces', faces, 'Vertices', vertices, 'FaceColor', 'interp', 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'FaceVertexCData', v); %,'EdgeAlpha', 1
axis equal
colormap(get_colormap(v)) %colormap(turbo) % colormap(flip(hot))     turbo flip(slanCM('guppy'))
sigma = std(v);
mu = mean(v);

if sigma > 0 && ignoreoutlier

    if ~isnumeric(ignoreoutlier)
        k = 3;
    else
        k = ignoreoutlier;
    end

    clim(gather([mu - k * sigma, mu + k * sigma]));
end

box off
axis off
end

function [left_faces, right_faces, left_vertices, right_vertices] = separate_hemispheres(faces, vertices)
% Number of vertices
num_vertices = size(vertices, 1);

% Step 1: Build adjacency matrix
adj_matrix = sparse(num_vertices, num_vertices);

for i = 1:size(faces, 1)
    v1 = faces(i, 1);
    v2 = faces(i, 2);
    v3 = faces(i, 3);

    % Add edges between vertices in the face
    adj_matrix(v1, v2) = 1;
    adj_matrix(v2, v1) = 1;
    adj_matrix(v2, v3) = 1;
    adj_matrix(v3, v2) = 1;
    adj_matrix(v1, v3) = 1;
    adj_matrix(v3, v1) = 1;
end

% Step 2: Find connected components
G = graph(adj_matrix);
components = conncomp(G);

% Assume we have two components, left and right hemispheres
left_component = components == 1;
right_component = components == 2;

% Step 3: Extract left and right vertices, along with their original indices
left_vertex_indices = find(left_component(:));
right_vertex_indices = find(right_component(:));

left_vertices = vertices(left_vertex_indices, :);
right_vertices = vertices(right_vertex_indices, :);
% idx_old=(1:length(vertices))';
    %{
    % Step 4: Create a mapping from old indices to new indices
    left_index_map = zeros(num_vertices, 1);
    right_index_map = zeros(num_vertices, 1);

    % Populate mapping with new indices
    left_index_map(left_vertex_indices) = 1:length(left_vertex_indices);
    right_index_map(right_vertex_indices) = 1:length(right_vertex_indices);

    % Step 5: Separate faces based on components and remap to new vertex indices
    left_face_mask = ismember(faces(:, 1), left_vertex_indices) & ...
                     ismember(faces(:, 2), left_vertex_indices) & ...
                     ismember(faces(:, 3), left_vertex_indices);
    right_face_mask = ismember(faces(:, 1), right_vertex_indices) & ...
                      ismember(faces(:, 2), right_vertex_indices) & ...
                      ismember(faces(:, 3), right_vertex_indices);

    left_faces = faces(left_face_mask, :);
    right_faces = faces(right_face_mask, :);

    % Remap old vertex indices in faces to the new indices
    left_faces = arrayfun(@(x) left_index_map(x), left_faces);
    right_faces = arrayfun(@(x) right_index_map(x), right_faces);
    %}
% left_vertex_indices
% right_vertex_indices
left_vertex_indices_new = 1:length(left_vertex_indices);
right_vertex_indices_new = 1:length(right_vertex_indices);
% faces
[~, loc] = ismember(faces, left_vertex_indices);
valid_entries = all(loc > 0, 2);
left_faces = left_vertex_indices_new(loc(valid_entries, :));

[~, loc] = ismember(faces, right_vertex_indices);
valid_entries = all(loc > 0, 2);
right_faces = right_vertex_indices_new(loc(valid_entries, :));

% idx_new=[left_vertex_indices;right_vertex_indices];
% faces=idx_new(faces);
% left_face_mask = ismember(faces(:, 1), left_vertex_indices) & ...
%                  ismember(faces(:, 2), left_vertex_indices) & ...
%                  ismember(faces(:, 3), left_vertex_indices);
% right_face_mask = ismember(faces(:, 1), right_vertex_indices) & ...
%                   ismember(faces(:, 2), right_vertex_indices) & ...
%                   ismember(faces(:, 3), right_vertex_indices);
% left_faces = faces(left_face_mask, :);
% right_faces = faces(right_face_mask, :);
end
