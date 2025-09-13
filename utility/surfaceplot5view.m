function [varargout] = surfaceplot5view(faces, vertices, vin, node, radius, isgpu, ignoreoutlier, opt)

if nargin < 7 || isempty(ignoreoutlier)
    ignoreoutlier = false;
end

if nargin < 6 || isempty(isgpu)
    isgpu = true;
end

% Define the radius of the smoothing kernel
if nargin < 5 || isempty(radius)
    radius = []; % replace with your desired radius
end

if nargin < 8 || isempty(opt)
    opt = struct();
end

opt = set_defaults(opt, 'isgpu', false);
opt = set_defaults(opt, 'isinflate', false);
opt = set_defaults(opt, 'inflated', false);

opt = set_defaults(opt, 'vfun', []);
opt = set_defaults(opt, 'isbinary', false);
opt = set_defaults(opt, 'reflective', false);
opt = set_defaults(opt, 'std_range', false);
opt = set_defaults(opt, 'tess', []);

if isempty(vertices) || isempty(faces)
    % Load cortex mesh: prefer opt.tess (struct with Vertices, Faces, Atlas).
    if ~isempty(opt.tess)
        tess = opt.tess;
    else
        this_file = mfilename('fullpath');
        repo_root = fileparts(fileparts(this_file)); % .../BiSCA
        mesh_path = fullfile(repo_root, 'data', 'sEEG', 'anatomy', 'anat', ...
            'tess_cortex_mid_low.mat');
        if exist(mesh_path, 'file') ~= 2
            error('surfaceplot5view:noMesh', ...
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

% split head vertices to left and right hemisphere also the faces of meshes
if ~iscell(vertices)
    %{
    % idx=vertices(:,2)>0;
    % vertices={vertices(idx,:),vertices(~idx,:)};
    % vin={vin(idx,:),vin(~idx,:)};
    % maxnode=max(faces,[],2);
    % minnode=min(faces,[],2);

    % idxleft=maxnode<=size(vertices{1},1);
    % idxright=minnode>size(vertices{1},1);
    % faces={faces(idxleft,:),faces(idxright,:)-size(vertices{1},1)+1};
    % faces={faces(all(faces <= size(leftHemisphereVertices, 1), 2), :),...
    %         faces(all(faces > size(leftHemisphereVertices, 1), 2), :)};
   % Identify vertices for each hemisphere
    midline=0;
    idx=vertices(:, 2) <= midline;
    leftHemisphereVertices = vertices(idx, :);
    rightHemisphereVertices = vertices(~idx, :);

    % Find indices of left and right hemisphere vertices
    leftVertexIndices  = ismember(vertices, leftHemisphereVertices, 'rows');
    rightVertexIndices = ismember(vertices, rightHemisphereVertices, 'rows');
    % Identify faces for each hemisphere
    leftHemisphereFaces = faces(all(leftVertexIndices(faces), 2), :);
    rightHemisphereFaces = faces(all(rightVertexIndices(faces), 2), :);

    vertices={leftHemisphereVertices,rightHemisphereVertices};
    vin={vin(idx,:),vin(~idx,:)};
    faces={leftHemisphereFaces,rightHemisphereFaces};
    %}
    face = faces; clear faces;
    vertice = vertices; clear vertices;
    [faces{1}, faces{2}, vertices{1}, vertices{2}] = separate_hemispheres(face, vertice);
    % vins{1}=vin(left_component);
    % vins{2}=vin(right_component);
    % vin=vins;clear vins vertice face;
end

if isgpu
    faces{1} = gpuArray(faces{1});
    faces{2} = gpuArray(faces{2});
    vertices{1} = gpuArray(vertices{1});
    vertices{2} = gpuArray(vertices{2});

    if iscell(vin)
        vin{1} = gpuArray(vin{1});
        vin{2} = gpuArray(vin{2});
    else
        vin = gpuArray(vin);
    end

    if nargin > 3 && ~isempty(node)

        if iscell(vin)
            node{1} = gpuArray(node{1});
            node{2} = gpuArray(node{2});
        else
            node = gpuArray(node);
        end

        radius = gpuArray(radius);
    end

end

if nargin > 3 && ~isempty(node)

    if iscell(vin)
        v{1} = node2vertice(vin{1}, vertices{1}, node{1}, radius);
        v{2} = node2vertice(vin{2}, vertices{2}, node{2}, radius);
    else
        v{1} = node2vertice(vin, vertices{1}, node, radius);
        v{2} = node2vertice(vin, vertices{2}, node, radius);
    end

else
    if iscell(vin)
        v{1} = vin{1};
        v{2} = vin{2};
    else
        % Support scalar per-vertex values passed as numeric vector.
        v{1} = vin;
        v{2} = vin;
    end
end


%{
% if nargout<2
% tiledlayout(2,2,"TileSpacing","compact")
% p.p1=surfpatch(faces{1},vertices{1},v{1});
% view(-90,0);
% p.p2=surfpatch(faces{1},vertices{1},v{1});
% view(90,0);
%
% p.p3=surfpatch(faces{1},vertices{1},v{1});
% view(0,90);
% p.p4=surfpatch(faces{2},vertices{2},v{2});
% view(0,90);
% c=colorbar;
% c.Location="south";
%
% p.p5=surfpatch(faces{2},vertices{2},v{2});
% view(90,0);
% p.p6=surfpatch(faces{2},vertices{2},v{2});
% view(-90,0);

% create the figure
% figure('Position',[1000,375,820,512]);
f=gcf;
%f.Position=[1000,375,820,512];
f.Position=[100	200	1520	960];
% define the positions of the subplots
left = 0.05;
bottom = 0.05;
hspacing = 0;
wspacing = 0;
width = (1-left*2-wspacing*2)/3;
height = (1-bottom*2-hspacing)/2;

% create the subplots
f1=subplot('Position', [left, bottom + height + hspacing, width, height]);
p.p1=surfpatch(faces{1},vertices{1},v{1},ignoreoutlier);view(-90,0);

f2=subplot('Position', [left , bottom , width, height]);
p.p2=surfpatch(faces{1},vertices{1},v{1},ignoreoutlier);view(90,0);

f3=subplot('Position', [left + 2*width + 2*wspacing, bottom + height + hspacing, width, height]);
p.p3=surfpatch(faces{2},vertices{2},v{2},ignoreoutlier);view(90,0);

f4=subplot('Position', [left + 2*width + 2*wspacing, bottom, width, height]);
p.p4=surfpatch(faces{2},vertices{2},v{2},ignoreoutlier);view(-90,0);

f5=subplot('Position', [left + width + wspacing, bottom+1/2*height, width, height]);
p.p5=surfpatch(faces{1},vertices{1},v{1},ignoreoutlier);view(0,90);hold on
p.p6=surfpatch(faces{2},vertices{2},v{2},ignoreoutlier);view(0,90);
c=colorbar;
c.Location="southoutside";

% twidth = 2*width*1/2;
% height = (1-bottom*2-hspacing)/2;

% set(f1,'position',[-(1/3)*width, -(1/4)*height+ height, 2*width, 2*height]);
% set(f2,'position',[-(1/3)*width, -(1/4)*height, 2*width, 2*height]);
% set(f3,'position',[-(1/3)*width+ 2*width, -(1/4)*height+height + hspacing, 2*width, 2*height]);
% set(f4,'position',[-(1/3)*width+ 2*width, -(1/4)*height, 2*width, 2*height]);
% set(f5,'position',[-(1/3)*width + 3/2*width, 2/3*height, width, height]);
% cpos=c.Position;
% c.Position=[cpos(1)+cpos(3)/4,cpos(2),cpos(3)/2,cpos(4)];
set(f1,'position',[-(59/100)*width, -(59/100)*height+ height, 2.5*width, 2.5*height]);
set(f2,'position',[-(59/100)*width, -(59/100)*height, 2.5*width, 2.5*height]);
set(f3,'position',[-(59/100)*width+ 2*width, -(59/100)*height+height + hspacing, 2.5*width, 2.5*height]);
set(f4,'position',[-(59/100)*width+ 2*width, -(59/100)*height, 2.5*width, 2.5*height]);
set(f5,'position',[-(59/100)*width + 1.59*width, 1/2*height, 1.3*width, 1.3*height]);
cpos=c.Position;
c.Position=[cpos(1)+cpos(3)/4,cpos(2),cpos(3)/2,cpos(4)];
%}

% end

if opt.inflated
    faces = tess_inflated.faces;
    vertices = tess_inflated.vertices;
end

clf
f = gcf;
f.Position = [100	200	500	320];
% define the positions of the subplots
left = 0.05;
bottom = 0.05;
hspacing = 0;
wspacing = 0;
width = (1 - left * 2 - wspacing * 2) / 3;
height = (1 - bottom * 2 - hspacing) / 2;
if isempty(node)
    idx_left = false(0, 1);
    all_x = [vertices{1}(:, 1); vertices{2}(:, 1)];
    all_y = [vertices{1}(:, 2); vertices{2}(:, 2)];
    all_z = [vertices{1}(:, 3); vertices{2}(:, 3)];
else
    idx_left = node(:, 2) > 0;
    all_x = [vertices{1}(:, 1); vertices{2}(:, 1); node(:, 1)];
    all_y = [vertices{1}(:, 2); vertices{2}(:, 2); node(:, 2)];
    all_z = [vertices{1}(:, 3); vertices{2}(:, 3); node(:, 3)];
end
xlimits = 1.2*[min(all_x), max(all_x)];
ylimits = 1.2*[min(all_y), max(all_y)];
zlimits = 1.2*[min(all_z), max(all_z)];


% zlim() requires increasing numeric limits
if ~isnumeric(zlimits) || numel(zlimits) ~= 2 || ~isfinite(zlimits(1)) || ~isfinite(zlimits(2)) || zlimits(2) <= zlimits(1)
    zlimits = 'auto';
end

factor = 1.8;
mksize = 20;

if isfield(opt, 'mksize') && ~isempty(opt.mksize)
    mksize = opt.mksize;
end

if ~isempty(opt.vfun)
    vin = opt.vfun(vin);
    v{1} = opt.vfun(v{1});
    v{2} = opt.vfun(v{2});
end

sigma = std(vin);
mu = mean(vin);
% climits = double([min(vin), max(vin)]);
if isfield(opt,'climits')&&~isempty(opt.climits)
    climits=opt.climits;
else
    climits = double([min([v{1};v{2}]), max([v{1};v{2}])]);
end
if ~isnumeric(climits) || numel(climits) ~= 2 || ~isfinite(climits(1)) || ~isfinite(climits(2)) || climits(2) <= climits(1)
    climits = 'auto';
end

% cmp=get_colormap(vin);
% cmp = 'hot';
% colormap(redbluecmap)
% cmp = redbluecmap;
cmp = plt.inferno;
if sigma > 0 && ~opt.isbinary && opt.std_range
    % colormap(turbo)% colormap(flip(hot))
    % clim(gather([mu-3*sigma,mu+3*sigma]));
    % climits = (gather([mu - 1 * sigma, mu + 1 * sigma]));
end

if opt.isbinary
    % cmp=jet(100);
    cmp = flip(autumn(100), 1);
    cmp = [cmp(1, :); cmp(end, :)];
end

if isfield(opt, 'cmp') && ~isempty(opt.cmp)
    cmp = opt.cmp;
end

% create the subplots
f1 = subplot('Position', [left, bottom + height + hspacing, width, height]);
% surfplot(faces{1}, vertices{1})
hSurf=surfpatch(faces{1}, vertices{1}, v{1}, ignoreoutlier); hold on;
% scatter3(node(idx_left, 1), node(idx_left, 2), node(idx_left, 3), mksize, vin(idx_left), 'filled', 'o');
% hSurf.Material = 'dull';

if opt.reflective

    if ischar(cmp) || isstring(cmp)
        cin = eval([cmp, '(', num2str(length(unique(vin))), ')']);
    elseif isnumeric(cmp) && size(cmp, 2) == 3
        cin = cmp;
    elseif iscell(cmp)
        cin = cmp;
    end
    [~,cidx]=map2color(vin,'NumColors',length(cin));
    cin = cin(cidx, :);
    reflectiveScatter3(node(idx_left, 1), node(idx_left, 2), node(idx_left, 3), mksize, cin);
else
    scatter3(node(idx_left, 1), node(idx_left, 2), node(idx_left, 3), mksize, vin(idx_left), 'filled', 'o');
end

view(180, 0);
xlim(factor * xlimits)
% ylim([0, ylimits(2)])
ylim([ylimits])
if ischar(zlimits)
    zlim('auto')
else
    zlim(zlimits)
end
if ischar(climits)
    clim('auto')
else
    clim(climits)
end
colormap(cmp)


f2 = subplot('Position', [left, bottom, width, height]);
% surfplot(faces{1}, vertices{1})
hSurf=surfpatch(faces{1}, vertices{1}, v{1}, ignoreoutlier); hold on;
% scatter3(node(idx_left, 1), node(idx_left, 2), node(idx_left, 3), mksize, vin(idx_left), 'filled', 'o');
% hSurf.Material = 'dull';

if opt.reflective

    if ischar(cmp) || isstring(cmp)
        cin = eval([cmp, '(', num2str(length(unique(vin))), ')']);
    elseif isnumeric(cmp) && size(cmp, 2) == 3
        cin = cmp;
    end
    [~,cidx]=map2color(vin,'NumColors',length(cin));
    cin = cin(cidx, :);
    reflectiveScatter3(node(idx_left, 1), node(idx_left, 2), node(idx_left, 3), mksize, cin);
else
    scatter3(node(idx_left, 1), node(idx_left, 2), node(idx_left, 3), mksize, vin(idx_left), 'filled', 'o');
end

view(0, 0);
xlim(factor * xlimits)
% ylim([0, ylimits(2)])
ylim([ylimits])
if ischar(zlimits)
    zlim('auto')
else
    zlim(zlimits)
end
if ischar(climits)
    clim('auto')
else
    clim(climits)
end
colormap(cmp)


f3 = subplot('Position', [left + 2 * width + 2 * wspacing, bottom + height + hspacing, width, height]);
% surfplot(faces{2}, vertices{2})
surfpatch(faces{2}, vertices{2}, v{2}, ignoreoutlier); hold on;
% scatter3(node(~idx_left, 1), node(~idx_left, 2), node(~idx_left, 3), mksize, vin(~idx_left), 'filled', 'o');
if opt.reflective

    if ischar(cmp) || isstring(cmp)
        cin = eval([cmp, '(', num2str(length(unique(vin))), ')']);
    elseif isnumeric(cmp) && size(cmp, 2) == 3
        cin = cmp;
    end
    [~,cidx]=map2color(vin,'NumColors',length(cin));
    cin = cin(cidx, :);
    reflectiveScatter3(node(~idx_left, 1), node(~idx_left, 2), node(~idx_left, 3), mksize, cin);
else
    scatter3(node(~idx_left, 1), node(~idx_left, 2), node(~idx_left, 3), mksize, vin(~idx_left), 'filled', 'o');
end

view(0, 0);
xlim(factor * (xlimits))
% ylim([ylimits(1), 0])
ylim([ylimits])
if ischar(zlimits)
    zlim('auto')
else
    zlim(zlimits)
end
if ischar(climits)
    clim('auto')
else
    clim(climits)
end
colormap(cmp)


f4 = subplot('Position', [left + 2 * width + 2 * wspacing, bottom, width, height]);
% surfplot(faces{2}, vertices{2})
hSurf=surfpatch(faces{2}, vertices{2}, v{2}, ignoreoutlier); hold on;
% scatter3(node(~idx_left, 1), node(~idx_left, 2), node(~idx_left, 3), mksize, vin(~idx_left), 'filled', 'o');
% hSurf.Material = 'dull';

if opt.reflective

    if ischar(cmp) || isstring(cmp)
        cin = eval([cmp, '(', num2str(length(unique(vin))), ')']);
    elseif isnumeric(cmp) && size(cmp, 2) == 3
        cin = cmp;
    end
    [~,cidx]=map2color(vin,'NumColors',length(cin));
    cin = cin(cidx, :);
    reflectiveScatter3(node(~idx_left, 1), node(~idx_left, 2), node(~idx_left, 3), mksize, cin);
else
    scatter3(node(~idx_left, 1), node(~idx_left, 2), node(~idx_left, 3), mksize, vin(~idx_left), 'filled', 'o');
end

view(-180, 0);
xlim(factor * (xlimits))
% ylim([ylimits(1), 0])
ylim([ylimits])
zlim(zlimits)
clim(climits)

f5 = subplot('Position', [left + width + wspacing, bottom+(1/2 * height), width, height]);
% surfplot(faces{1}, vertices{1})
% surfplot(faces{2}, vertices{2})
hSurf=surfpatch(faces{1}, vertices{1}, v{1}, ignoreoutlier); hold on;
% hSurf.Material = 'dull';

hSurf=surfpatch(faces{2}, vertices{2}, v{2}, ignoreoutlier);
% hSurf.Material = 'dull';

% scatter3(node(:, 1), node(:, 2), node(:, 3), mksize, vin, 'filled', 'o');
if opt.reflective

    if ischar(cmp) || isstring(cmp)
        cin = eval([cmp, '(', num2str(length(unique(vin))), ')']);
    elseif isnumeric(cmp) && size(cmp, 2) == 3
        cin = cmp;
    end
    [~,cidx]=map2color(vin,'NumColors',length(cin));

    cin = cin(cidx, :);
    reflectiveScatter3(node(:, 1), node(:, 2), node(:, 3), mksize, cin);
else
    scatter3(node(:, 1), node(:, 2), node(:, 3), mksize, vin, 'filled', 'o');
end

view(-90, 90);
c = colorbar;
c.Location = "southoutside";
c.FontSize=14;
xlim(xlimits)
ylim(ylimits)
zlim(zlimits)
clim(climits)
colormap(cmp)
% c.Location = "manual";
% c.Location=[0.396166666666667,0.145094086021505,0.195,0.0585];
 c.Position=[0.355,0.15,0.3,0.05];

set(f1, 'position', [- (33/100) * width, - (35/100) * height + height, 2.0 * width, 2.0 * height]);
set(f2, 'position', [- (33/100) * width, - (35/100) * height, 2 * width, 2 * height]);
set(f3, 'position', [- (33/100) * width + 2 * width, - (35/100) * height + height + hspacing, 2.0 * width, 2.0 * height]);
set(f4, 'position', [- (33/100) * width + 2 * width, - (35/100) * height, 2.0 * width, 2.0 * height]);
set(f5, 'position', [- (36/100) * width + 1.4 * width, 1/2 * height, 1.3 * width, 1.3 * height]);
cpos = c.Position;
c.Position = [cpos(1) + cpos(3) / 4, cpos(2), cpos(3) / 2, cpos(4)];

if nargout == 0
    varargout = {};
end

if nargout == 1
    varargout{1} = p;
else
    varargout{1} = [];
end

if nargout >= 2

    if size(vin, 2) == 1
        varargout{2} = v;
    else
        varargout{2} = [];
    end

end

end

function p = surfpatch(faces, vertices, v, ignoreoutlier)
p = patch('Faces', faces, 'Vertices', vertices, 'FaceColor', 'interp', 'FaceAlpha', 0.7, 'EdgeColor', 'none', 'FaceVertexCData', v); %,'EdgeAlpha', 1
axis equal
% colormap(get_colormap(v)) %colormap(turbo) % colormap(flip(hot))     turbo flip(slanCM('guppy'))
colormap(plt.inferno)

sigma = std(v);
mu = mean(v);

% if sigma > 0 && ignoreoutlier
% 
%     if ~isnumeric(ignoreoutlier)
%         k = 3;
%     else
%         k = ignoreoutlier;
%     end
% 
%     clim(gather([mu - k * sigma, mu + k * sigma]));
% end

box off
axis off
end

% figure(1);clf;hold on
% if nargout<2
%
%
%     % v(v<mu-3*sigma)=mu-3*sigma;
%     % v(v>mu+3*sigma)=mu+3*sigma;
%

%
%     if size(vin,2)>1
%         for i=1:size(vin,2)
%             if size(vin,1)~=size(vertices)
%                 v=weights*vin(i,:);
%             else
%                 v=vin(i,:);
%             end
%             set(p,'FaceVertexCData',v);
%             sigma=std(v);
%             mu=mean(v);
%             clim([mu-3*sigma,mu+3*sigma]);
%             pause(0.2)
%         end
%     end
%
%
% end

function v = node2vertice(vin, vertices, node, radius)
% v=NwSmooth(node,vin,radius,vertices);
if nargin < 4 || isempty(radius)
    radius = max(max(vertices, [], 1) - min(vertices, [], 1) / 2);
    % radius=radius*ones(size(vertices,2),1);
    m = 40;
    hlist = get_hList(m, [radius / 100, radius * 5], @logspace);
    [~, nw_gcv] = NwSmoothGCV(gather(node), gather(vin), gather(hlist));
    % [nw_gcv]=cvNwSmooth(node,vin,node,m,hlist);
    v = NwSmooth(node, vin, nw_gcv.hmin, vertices);
else
    v = NwSmooth(node, vin, radius, vertices);
end

end

function reflectiveScatter3(x, y, z, mksize, colors)
% Custom function to plot reflective points as small spheres.
% Inputs:
% x, y, z - coordinates of the points
% mksize - size of the markers
% colors - color matrix (Nx3) for each point

% Check if colors are provided; if not, default to a single color
if nargin < 5 || isempty(colors)
    colors = repmat([0.5, 0.5, 0.5], length(x), 1); % Default gray color
end

% Create a sphere with 10 subdivisions for smoother visualization
[sx, sy, sz] = sphere(10);
mksize = mksize ./ 5000;
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
% camlight headlight; % Add a camera light for better visualization

% Final touch to make the scene look better
axis equal;
hold off;
end

% function v = node2vertice(vin, vertices, node, radius)
% if size(vin,1)~=size(vertices)
%     % Compute the distance between the ChannelPosition and the Nodes
%     distances = pdist2(vertices,node);
%
%     % Compute the weights of the smoothing kernel
%     weights = exp(-(distances.^2)/(2*radius^2));
%
%     % Normalize the weights
%     weights = weights./sum(weights,2);
%
%     % Compute the smoothed value at the ChannelPosition
%     v = weights*vin(:,1);
%     if isempty(find(isnan(vin(:,1)),1))&&~isempty(find(isnan(v),1))
%         warning('has missing value, enlarge your radius or use leadfield to smooth')
%     end
% else
%     v=vin;
% end
% end

%{
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
    left_vertex_indices = find(left_component);
    right_vertex_indices = find(right_component);

    left_vertices = vertices(left_vertex_indices, :);
    right_vertices = vertices(right_vertex_indices, :);

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
end
%}
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

%{

function [left_faces, right_faces, left_vertices, right_vertices, left_component, right_component] = separate_hemispheres(faces, vertices)
% Number of vertices
num_vertices = size(vertices, 1);

% Step 1: Build adjacency matrix
adj_matrix = sparse(num_vertices, num_vertices);

for i = size(faces, 1):-1:1
    v1 = faces(i, 1);
    v2 = faces(i, 2);
    v3 = faces(i, 3);

    % Add edges between vertices in the face
    adj_matrix(v1, v2) = 1; %#ok<*SPRIX>
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

% Step 3: Separate faces based on components
left_faces = faces(ismember(faces(:,1), find(left_component)) & ...
    ismember(faces(:,2), find(left_component)) & ...
    ismember(faces(:,3), find(left_component)), :);

right_faces = faces(ismember(faces(:,1), find(right_component)) & ...
    ismember(faces(:,2), find(right_component)) & ...
    ismember(faces(:,3), find(right_component)), :);

% Step 4: Separate vertices based on components
left_vertices = vertices(left_component, :);
right_vertices = vertices(right_component, :);
end

%}
