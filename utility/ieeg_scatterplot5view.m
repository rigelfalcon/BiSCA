function ieeg_scatterplot5view(iEEG, varname, opt)
% IEEG_SCATTERPLOT5VIEW  5-view scatter plot of iEEG channel values on cortex
%   ieeg_scatterplot5view(iEEG, varname, opt)
%
%   Inputs:
%     iEEG    - struct with fields: subjinfo, Nsubj, cortex.tess
%     varname - field name in subjinfo(i) containing per-channel values
%     opt     - options struct (optional):
%               .isbinary    - binary colormap mode (default: false)
%               .reflective  - use reflective spheres (default: false)
%               .vfun        - function handle applied to values (default: [])
%               .cmp         - custom colormap (default: [])
%               .mksize      - marker size (default: 30)

if nargin < 3 || isempty(opt)
    opt = struct();
end
opt = set_defaults(opt, 'isbinary', false);
opt = set_defaults(opt, 'reflective', false);
opt = set_defaults(opt, 'vfun', []);
opt = set_defaults(opt, 'cmp', []);
opt = set_defaults(opt, 'mksize', 30);

%% Load cortex mesh from iEEG struct
tess = iEEG.cortex.tess;
vertices{1} = tess.Vertices(tess.Atlas(2).Scouts(1).Vertices, :);
vertices{2} = tess.Vertices(tess.Atlas(2).Scouts(2).Vertices, :);
faces{1} = tess.Faces(tess.Faces(:,1) < tess.Atlas(2).Scouts(2).Seed, :);
faces{2} = tess.Faces(tess.Faces(:,1) >= tess.Atlas(2).Scouts(2).Seed, :);
faces{2} = faces{2} - tess.Atlas(2).Scouts(2).Seed + 1;

%% Collect channel positions and values
ichall = 0;
for isubj = iEEG.Nsubj:-1:1
    Nch = iEEG.subjinfo(isubj).NumChannel;
    for ich = Nch:-1:1
        ichall = ichall + 1;
        node(ichall, :) = iEEG.subjinfo(isubj).ChannelPositionMidBst(ich, :);
        vin(ichall, :) = iEEG.subjinfo(isubj).(varname)(ich);
        if isempty(iEEG.subjinfo(isubj).(varname)(ich))
            error('lost value');
        end
    end
end

%% Apply value transform
if ~isempty(opt.vfun)
    vin = opt.vfun(vin);
end

%% Compute color limits
mksize = opt.mksize;
sigma = std(vin);
mu = mean(vin);
climits = double([min(vin), max(vin)]);
cmp = [];
if sigma > 0 && ~opt.isbinary
    climits = gather([mu - 1*sigma, mu + 1*sigma]);
end
if opt.isbinary
    cmp = flip(autumn(100), 1);
    cmp = [cmp(1,:); cmp(end,:)];
end
if isfield(opt, 'cmp') && ~isempty(opt.cmp) && ~opt.isbinary
    cmp = opt.cmp;
end

%% Setup figure and layout
clf
f = gcf;
Position = [200, 200, 660, 410];
f.Position = Position;

left = 0.05;
bottom = 0.05;
hspacing = 0;
wspacing = 0;
width = (1 - left*2 - wspacing*2) / 3;
height = (1 - bottom*2 - hspacing) / 2;
idx_left = node(:,2) > 0;
xlimits = 1.2 * [min(min(tess.Vertices(:,1), min(node(:,1)))), max(max(tess.Vertices(:,1), max(node(:,1))))];
ylimits = 1.2 * [min(min(tess.Vertices(:,2), min(node(:,2)))), max(max(tess.Vertices(:,2), max(node(:,2))))];
zlimits = 1.2 * [min(min(tess.Vertices(:,3), min(node(:,3)))), max(max(tess.Vertices(:,3), max(node(:,3))))];
factor = 1.8;

%% Helper: get color index for reflective mode
    function cin = get_cin(cmp_in, vin_in)
        if ischar(cmp_in) || isstring(cmp_in)
            cin = eval([cmp_in, '(', num2str(length(unique(vin_in))), ')']);
            cin = cin(vin_in, :);
        elseif isnumeric(cmp_in) && size(cmp_in, 2) == 3
            cin = cmp_in(vin_in, :);
        else
            cmp_fb = plt.viridis(256);
            vn = (vin_in - min(vin_in)) / max(range(vin_in), eps);
            cin = cmp_fb(max(1, min(256, round(vn * 255) + 1)), :);
        end
    end

%% Pre-compute colors for reflective mode (once for all subplots)
if opt.reflective
    cin_all = get_cin(cmp, vin);
end

%% Subplot 1: Left hemisphere, posterior view
f1 = subplot('Position', [left, bottom + height + hspacing, width, height]);
surfplot(faces{1}, vertices{1});
if opt.reflective
    reflectiveScatter3(node(idx_left,1), node(idx_left,2), node(idx_left,3), mksize, cin_all(idx_left,:));
else
    scatter3(node(idx_left,1), node(idx_left,2), node(idx_left,3), mksize, vin(idx_left), 'filled', 'o');
end
xlim(factor*xlimits); ylim(ylimits); zlim(zlimits); clim(climits);
view(180, 0); camlight headlight; camlight headlight; colormap(cmp);

%% Subplot 2: Left hemisphere, anterior view
f2 = subplot('Position', [left, bottom, width, height]);
surfplot(faces{1}, vertices{1});
if opt.reflective
    reflectiveScatter3(node(idx_left,1), node(idx_left,2), node(idx_left,3), mksize, cin_all(idx_left,:));
else
    scatter3(node(idx_left,1), node(idx_left,2), node(idx_left,3), mksize, vin(idx_left), 'filled', 'o');
end
xlim(factor*xlimits); ylim(ylimits); zlim(zlimits); clim(climits);
view(0, 0); camlight headlight; camlight headlight; colormap(cmp);

%% Subplot 3: Right hemisphere, anterior view
f3 = subplot('Position', [left + 2*width + 2*wspacing, bottom + height + hspacing, width, height]);
surfplot(faces{2}, vertices{2});
if opt.reflective
    reflectiveScatter3(node(~idx_left,1), node(~idx_left,2), node(~idx_left,3), mksize, cin_all(~idx_left,:));
else
    scatter3(node(~idx_left,1), node(~idx_left,2), node(~idx_left,3), mksize, vin(~idx_left), 'filled', 'o');
end
xlim(factor*xlimits); ylim(ylimits); zlim(zlimits); clim(climits);
view(0, 0); camlight headlight; camlight headlight; colormap(cmp);

%% Subplot 4: Right hemisphere, posterior view
f4 = subplot('Position', [left + 2*width + 2*wspacing, bottom, width, height]);
surfplot(faces{2}, vertices{2});
if opt.reflective
    reflectiveScatter3(node(~idx_left,1), node(~idx_left,2), node(~idx_left,3), mksize, cin_all(~idx_left,:));
else
    scatter3(node(~idx_left,1), node(~idx_left,2), node(~idx_left,3), mksize, vin(~idx_left), 'filled', 'o');
end
xlim(factor*xlimits); ylim(ylimits); zlim(zlimits); clim(climits);
view(-180, 0); camlight headlight; camlight headlight; colormap(cmp);

%% Subplot 5: Top-down view (both hemispheres)
f5 = subplot('Position', [left + width + wspacing, bottom + 1/2*height, width, height]);
surfplot(faces{1}, vertices{1});
surfplot(faces{2}, vertices{2});
if opt.reflective
    reflectiveScatter3(node(:,1), node(:,2), node(:,3), mksize, cin_all);
else
    scatter3(node(:,1), node(:,2), node(:,3), mksize, vin, 'filled', 'o');
end
xlim(xlimits); ylim(ylimits); zlim(zlimits); clim(climits);
view(-90, 90); camlight headlight; camlight headlight; colormap(cmp);

cb = colorbar;
cb.Location = 'southoutside';
cb.TickLabelsMode = 'manual';
cb.FontSize = 8;
cb.TickLabels = {'L+G', 'L+NG', 'NL+G', 'NL+NG'};
cb.Ticks = [cb.Limits(1) + (cb.Limits(2)-cb.Limits(1))/8 : (cb.Limits(2)-cb.Limits(1))/4 : cb.Limits(2)+cb.Limits(1)];
cb.Position = [0.38, 0.16, 0.26, 0.018];

%% Adjust subplot positions
set(f1, 'position', [-(33/100)*width, -(35/100)*height + height, 2.0*width, 2.0*height]);
set(f2, 'position', [-(33/100)*width, -(35/100)*height, 2*width, 2*height]);
set(f3, 'position', [-(33/100)*width + 2*width, -(35/100)*height + height + hspacing, 2.0*width, 2.0*height]);
set(f4, 'position', [-(33/100)*width + 2*width, -(35/100)*height, 2.0*width, 2.0*height]);
set(f5, 'position', [-(36/100)*width + 1.4*width, 1/2*height, 1.3*width, 1.3*height]);

end


function surfplot(faces, vertices)
% Render semi-transparent cortex surface
hold on; box off; axis off; axis equal;
patch('Faces', faces, 'Vertices', vertices, ...
    'FaceColor', 'flat', 'FaceAlpha', 0.1, ...
    'EdgeColor', 0.1*ones(1,3), ...
    'FaceVertexCData', 0.5*ones(length(vertices), 3));
shading interp;
material metal;
lighting phong;
end


function reflectiveScatter3(x, y, z, mksize, colors)
% Plot reflective points as small spheres
if nargin < 5 || isempty(colors)
    colors = repmat([0.5, 0.5, 0.5], length(x), 1);
end
[sx, sy, sz] = sphere(10);
mksize = mksize ./ 13000;
hold on;
for i = 1:length(x)
    surf(x(i) + sx * mksize, y(i) + sy * mksize, z(i) + sz * mksize, ...
        'FaceColor', colors(i, :), 'EdgeColor', 'none');
end
material metal;
lighting phong;
axis equal;
hold off;
end
