%% fig3_ieeg_spatial_freq_statistic
% Script:   fig3_ieeg_peak_frequency.m
% Figure:   Fig. 3e-f (iEEG peak frequency)
% Purpose:  Generate spatial and frequency summary statistics for iEEG.
% Usage:    matlab -batch "run('./setup.m'); run('./example/BiSCA/fig3_ieeg_peak_frequency.m')"
%

close all
clearvars -except stat_bh

% Get repo root
this_file = mfilename('fullpath');
candidate_roots = {
    fileparts(fileparts(fileparts(this_file)));  % 3-deep
    fileparts(fileparts(this_file));             % 2-deep
};
repo_root = [];
for i = 1:numel(candidate_roots)
    if exist(fullfile(candidate_roots{i}, 'setup.m'), 'file') == 2
        repo_root = candidate_roots{i};
        break;
    end
end
assert(~isempty(repo_root), 'Cannot find repository root containing setup.m');

% Load data if not already loaded
if ~exist('iEEG', 'var') || ~isfield(iEEG.subjinfo(1), 'bixialpha')
    filepath = fullfile(repo_root, 'data', 'result', 'iEEG_bixialpha.mat');
    load(filepath);
end

if ~exist('stat_bh', 'var')
    filepath = fullfile(repo_root, 'data', 'result', 'stat_bh_summary.mat');
    load(filepath);
end

stat_all = stat_bh;
alpha = stat_all.raw.median(1).p_threshold;

% Extract channel and bixialpha data
Nch_subj = cat_struct_fields(iEEG.subjinfo, 'NumChannel', false, true);
idx_subj_ch = gen_idx_subj_ch(Nch_subj);
bixialpha_all = cat_struct_fields(iEEG.subjinfo, 'bixialpha', false, true);
%%
if ~exist('stat_all', 'var')
    stat_all.max_abs_bc_f = [];
    stat_all.max_abs_bc = [];
    % Initialize stat_all structure
    for idx = length(idx_subj_ch):-1:1
        bixialpha = bixialpha_all(idx);
        hos.bc = bixialpha.bc;
        hos.f  = bixialpha.f;
        hos.Ns = bixialpha.Ns;

        % Extract bc matrix for the current channel
        bc_matrix = hos.bc;

        %{
        % Calculate maximum absolute value and its location
        [max_val, linear_idx] = max(abs(bc_matrix(:)));
        [i_idx, j_idx] = ind2sub(size(bc_matrix), linear_idx);

        % Store results
        stat_all.max_abs_bc(idx, :) = max_val;
        stat_all.median_abs_bc(idx, :) = median(abs(bc_matrix(:)), "all", "omitmissing");
        stat_all.max_abs_bc_f(idx, :) = [bixialpha.f(i_idx), bixialpha.f(j_idx)];
        %}
        bc_dig = abs(full2diag(bc_matrix));
        [max_abs_bc_diag, idx_max] = max(bc_dig, [], 1);
        frange = linspace(0, 50, size(bc_dig, 1));
        f = frange(idx_max);

        stat_all.max_abs_bc_f = [f(:); stat_all.max_abs_bc_f];
        stat_all.max_abs_bc = [max_abs_bc_diag(:); stat_all.max_abs_bc];
    end

end

%% Calculate probabilities for GL, NGL, GNL, NGNL
% Get HG and HL values
H_G = [stat_all.raw.median.h];
H_L = [stat_all.raw.max.h];

% Calculate combinations
H_GL = ~H_G & ~H_L;
H_NGL = H_G & ~H_L;
H_GNL = ~H_G & H_L;
H_NGNL = H_G & H_L;

% Calculate probabilities
p_GL = sum(H_GL) / length(H_GL);
p_NGL = sum(H_NGL) / length(H_GL);
p_GNL = sum(H_GNL) / length(H_GL);
p_NGNL = sum(H_NGNL) / length(H_GL);

disp(['p_GL: ', num2str(p_GL), ' p_NGL: ', num2str(p_NGL), ' p_GNL: ', num2str(p_GNL), ' p_NGNL: ', num2str(p_NGNL)]);

%% Plot frequency probabilities
% Get frequency data
f = linspace(0, 50, size(stat_all.raw.median(1).bcth, 1));
[fx, fy] = meshgrid(f, f);

bcth_NG = cellfun(@(x) full(x), {stat_all.raw.median.bcth}, 'UniformOutput', false)';
N = size(bcth_NG{1}, 1);
bcth_NG = cat(3, bcth_NG{:});
bcth_NL = cellfun(@(x) full(x), {stat_all.raw.max.bcth}, 'UniformOutput', false)';
bcth_NL = cat(3, bcth_NL{:});

% element-wise
h_GL = logical(~abs(bcth_NG) & ~abs(bcth_NL));
h_NGL = logical(abs(bcth_NG) & ~abs(bcth_NL));
h_GNL = logical(~abs(bcth_NG) & abs(bcth_NL));
h_NGNL = logical(abs(bcth_NG) & abs(bcth_NL));
h_NG = logical(abs(bcth_NG));
h_NL = logical(abs(bcth_NL));

% Normalize probabilities
p_freq_GL = sum(h_GL, 3) ./ size(h_GL, 3);
p_freq_NGL = sum(h_NGL, 3) ./ size(h_NGL, 3);
p_freq_GNL = sum(h_GNL, 3) ./ size(h_GNL, 3);
p_freq_NGNL = sum(h_NGNL, 3) ./ size(h_NGNL, 3);
p_freq_NG = sum(h_NG, 3) ./ size(h_NG, 3);
p_freq_NL = sum(h_NL, 3) ./ size(h_NL, 3);

% Apply mask for valid frequencies
idx = (fx + fy) < f(end);
p_freq_GL(~idx) = NaN;
p_freq_NGL(~idx) = NaN;
p_freq_GNL(~idx) = NaN;
p_freq_NGNL(~idx) = NaN;
p_freq_NG(~idx) = NaN;
p_freq_NL(~idx) = NaN;
p_freq_GL(1, :) = NaN;
p_freq_NGL(1, :) = NaN;
p_freq_GNL(1, :) = NaN;
p_freq_NGNL(1, :) = NaN;
p_freq_NG(1, :) = NaN;
p_freq_NL(1, :) = NaN;
p_freq_GL(:, 1) = NaN;
p_freq_NGL(:, 1) = NaN;
p_freq_GNL(:, 1) = NaN;
p_freq_NGNL(:, 1) = NaN;
p_freq_NG(:, 1) = NaN;
p_freq_NL(:, 1) = NaN;

%% Plot frequency probabilities
figure(401)
set(gcf, 'Position', [197.6667, 991, 450, 350]);

% Unified linear color limits across all 4 bifrequency panels (iEEG internal)
linear_data_max = max([p_freq_GL(:); p_freq_NGL(:); p_freq_GNL(:); p_freq_NGNL(:)], [], 'omitmissing');
fprintf('[fig3_ieeg] linear data max = %.4f\n', linear_data_max);
clims = [0, linear_data_max];
tiledlayout(2, 2) %,'TileSpacing', 'compact'
nexttile
surfbc(fx, fy, p_freq_GL)
clim(clims)
view(2)
colormap(plt.inferno(256))
colorbar

nexttile
surfbc(fx, fy, p_freq_NGL)
clim(clims)
view(2)
colormap(plt.inferno(256))
colorbar

nexttile
surfbc(fx, fy, p_freq_GNL)
clim(clims)
view(2)
colormap(plt.inferno(256))
colorbar

nexttile
surfbc(fx, fy, p_freq_NGNL)
clim(clims)
view(2)
colormap(plt.inferno(256))
colorbar

% Save frequency figure
filename = fullfile(repo_root, 'figure', 'fig3', 'fig3_iEEG_freq_stat.png');
test_folder(filename);
exportgraphics(gcf, filename, 'Resolution', 600);
if ispc; print(gcf, strrep(filename, '.png', '.emf'), '-dmeta'); end


%%
% Plot scatter surface


%%
% Plot scatter surface
opt.thfactor = 3; % Outlier removal threshold
opt.vfun = []; % No custom value function for scatter plots
opt.inflated = false;
opt.reflective = false;
opt.smooth_factor = 3;
FontSize = 12;

figure(405)
set(gcf, 'Position', [100, 99.6667, 700, 565]);
clf;

% Create tiledlayout with padding to leave space for colorbar
t = tiledlayout(2, 2, 'TileSpacing', 'none', 'Padding', 'loose');

% Unified color scale [0,1] for surfscatter cross-panel comparability
clims_surf = [0, 1];

ax1 = nexttile;
V = H_GL';
plot_scattersurf(iEEG, V, opt);
clim(clims_surf)
c1 = colorbar(ax1, 'Orientation', 'vertical');
c1.FontSize = FontSize;

% Second subplot
ax2 = nexttile;
V = H_NGL';
plot_scattersurf(iEEG, V, opt);
clim(clims_surf)
c2 = colorbar(ax2, 'Orientation', 'vertical');
c2.FontSize = FontSize;

% Third subplot
ax3 = nexttile;
V = H_GNL';
plot_scattersurf(iEEG, V, opt);
clim(clims_surf)
c3 = colorbar(ax3, 'Orientation', 'vertical');
c3.FontSize = FontSize;

% Fourth subplot
ax4 = nexttile;
V = H_NGNL';
plot_scattersurf(iEEG, V, opt);
clim(clims_surf)
c4 = colorbar(ax4, 'Orientation', 'vertical');
c4.FontSize = FontSize;

% Shorten colorbars to 60% of axis height
axes_list = [ax1, ax2, ax3, ax4];
cbars = [c1, c2, c3, c4];
fraction = 0.6; % Height fraction (60% of axis height)

for i = 1:4
    ax = axes_list(i);
    cbar = cbars(i);

    % Get current positions (normalized units)
    axPos = ax.Position;
    cPos = cbar.Position;

    % Calculate new height and vertical position
    newHeight = axPos(4) * fraction;
    newBottom = axPos(2) + (axPos(4) - newHeight) / 2;

    % Adjust colorbar position (keep width/horizontal position)
    cbar.Position = [cPos(1), newBottom, cPos(3), newHeight];
end

% Save figure
filename = fullfile(repo_root, 'figure', 'fig3', 'fig3_iEEG_surfscatter.png');
test_folder(filename);
exportgraphics(gcf, filename, 'Resolution', 600);
if ispc; print(gcf, strrep(filename, '.png', '.emf'), '-dmeta'); end

disp('Figure 3 generation: proceeding to EEG topography.');

%% EEG topography plots (uses local copied HeadModel + Channel)
hm_path = fullfile(repo_root, 'data', 'sEEG', 'anatomy', 'data', 'MNI19ch_Compare_Source', 'ScalpEEG', 'headmodel_surf_openmeeg.mat');
ch_path = fullfile(repo_root, 'data', 'sEEG', 'anatomy', 'data', 'MNI19ch_Compare_Source', 'ScalpEEG', 'channel_10-20_19.mat');
if ~exist(hm_path, 'file') || ~exist(ch_path, 'file')
    warning('Anatomy headmodel or channel file not found. Skipping scalp projection section.');
else
HeadModel = load(hm_path);
Channel = load(ch_path);
Lvj = bst_gain_orient(HeadModel.Gain, HeadModel.GridOrient);
%%
J_H_GL = electrodes2surface(iEEG, H_GL', opt);
V_H_GL = Lvj * J_H_GL;

J_H_NGL = electrodes2surface(iEEG, H_NGL', opt);
V_H_NGL = Lvj * J_H_NGL;

J_H_GNL = electrodes2surface(iEEG, H_GNL', opt);
V_H_GNL = Lvj * J_H_GNL;

J_H_NGNL = electrodes2surface(iEEG, H_NGNL', opt);
V_H_NGNL = Lvj * J_H_NGNL;

figure(406)
set(gcf, "Position", [458.3333333333333, 476.3333333333333, 400, 350]);
clims = [min([V_H_GL(:); V_H_NGL(:); V_H_GNL(:); V_H_NGNL(:)]) max([V_H_GL(:); V_H_NGL(:); V_H_GNL(:); V_H_NGNL(:)])];
tiledlayout(2, 2)
nexttile

plot_topography(V_H_GL, 'ChannelLabel', {Channel.Channel.Name});
clim(clims)
colorbar

nexttile

plot_topography(V_H_NGL, 'ChannelLabel', {Channel.Channel.Name});
clim(clims)
colorbar

nexttile

plot_topography(V_H_GNL, 'ChannelLabel', {Channel.Channel.Name});
clim(clims)
colorbar

nexttile

plot_topography(V_H_NGNL, 'ChannelLabel', {Channel.Channel.Name});
clim(clims)
colorbar

end % if anatomy files exist

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_scattersurf(iEEG, V, opt)
opt.reflective = false;
ChannelPos = cat_struct_fields(iEEG.subjinfo, 'ChannelPositionMidBst', false, true);

idx_bad = isnan(V) | isinf(V);
ChannelPos = ChannelPos(~idx_bad, :);
V = V(~idx_bad);

opt.inflated = false;
opt.tess = iEEG.cortex.tess;
opt.radius = max(max(iEEG.cortex.tess.Vertices, [], 1) - min(iEEG.cortex.tess.Vertices, [], 1)) * opt.smooth_factor;

surfscatter([], [], V, ChannelPos, opt);

colormap(plt.inferno(256))
end

%%
function J = electrodes2surface(iEEG, V, opt)
opt.inflated = false;
opt.radius = max(max(iEEG.cortex.tess.Vertices, [], 1) - min(iEEG.cortex.tess.Vertices, [], 1)) * 2;
ChannelPos = cat_struct_fields(iEEG.subjinfo, 'ChannelPositionMidBst', false, true);

idx_bad = isnan(V) | isinf(V);
ChannelPos = ChannelPos(~idx_bad, :);
V = V(~idx_bad);

J = node2vertice(V, iEEG.cortex.tess.Vertices, ChannelPos, opt.radius);
end

function plot_topography(channel_data, varargin)
% PLOT_TOPGRAPHY Create EEG channel topography plot with log color scaling
%   plot_topography(channel_data, 'Param', value, ...)
%
%   Inputs:
%   channel_data - Vector of values per channel (natural scale values)
%
%   Optional Parameters:
%   'Colormap'    - Custom colormap (default: viridis)
%   'CLim'        - Color axis limits in natural scale (default: [min max] of data)
%   'LogScale'    - Use log color scale (default: true)
%   'SigThreshold'- Significance threshold to mark on colorbar (natural scale)
%   'SavePath'    - Directory to save figure
%   'FileName'    - Base filename for saving
%   'Visible'     - Figure visibility ('on' or 'off')
%   'ChannelLabel'- Cell array of channel labels to include
%   'CBTitle'     - Colorbar title (default: 'Value')
%%

if ~exist('topoplot', 'file')
    eeglab('nogui')
end

%% Parameter parsing
p = inputParser;
addRequired(p, 'channel_data', @isnumeric);
addParameter(p, 'Colormap', plt.inferno(256), @(x) isnumeric(x) && size(x, 2) == 3);
addParameter(p, 'CLim', [], @isnumeric);
addParameter(p, 'LogScale', false, @islogical);
addParameter(p, 'SigThreshold', [], @isnumeric);
addParameter(p, 'SavePath', pwd, @ischar);
addParameter(p, 'FileName', 'topography_plot', @ischar);
addParameter(p, 'Visible', 'on', @(x) ismember(lower(x), {'on', 'off'}));
addParameter(p, 'ChannelLabel', [], @iscell);
addParameter(p, 'CBTitle', 'T-value', @ischar);
parse(p, channel_data, varargin{:});

%% Channel configuration
chanlocs = load_standard_channels();

if ~isempty(p.Results.ChannelLabel)
    idx = find_char({chanlocs.labels}, p.Results.ChannelLabel);
    chanlocs = chanlocs(idx(idx ~= 0));
end

% Validate input dimensions
if numel(channel_data) ~= numel(chanlocs)
    error('Input data (%d channels) doesn''t match configuration (%d channels)', ...
        numel(channel_data), numel(chanlocs));
end

%% Prepare data for plotting
if p.Results.LogScale
    plot_data = log10(channel_data);

    if ~isempty(p.Results.SigThreshold)
        sig_thresh_log = log10(p.Results.SigThreshold);
    end

else
    plot_data = channel_data;

    if ~isempty(p.Results.SigThreshold)
        sig_thresh_log = p.Results.SigThreshold;
    end

end

% Set color limits
if isempty(p.Results.CLim)
    cmin = min(channel_data(:));
    cmax = max(channel_data(:));
else
    cmin = p.Results.CLim(1);
    cmax = p.Results.CLim(2);
end

if p.Results.LogScale
    clims = log10([cmin cmax]);
else
    clims = [cmin cmax];
end

%% Generate topography plot
topoplot(plot_data, chanlocs, ...
    'maplimits', clims, ...
    'colormap', p.Results.Colormap, ...
    'electrodes', 'on', ...
    'shading', 'interp', ...
    'whitebk', 'on');
clim(clims)

% Handle log scale tick labels
if p.Results.LogScale
    % Get current tick marks in log space
    ticks = get(cb, 'Ticks');

    % Convert to natural scale for labels
    tick_labels = arrayfun(@(x) sprintf('%.2f', 10 .^ x), ticks, 'UniformOutput', false);
    set(cb, 'TickLabels', tick_labels);
end

end

function v = node2vertice(vin, vertices, node, radius)
if nargin < 4 || isempty(radius)
    radius = max((max(vertices, [], 1) - min(vertices, [], 1)) / 2);
    m = 100;
    hlist = get_hList(m, [radius / 100, radius * 5], @logspace); %logspace
    [~, nw_gcv] = NwSmoothGCVOriginalScaleLogGaussian(gather(node), gather(vin), gather(hlist));
    v = NwSmooth(node, vin, nw_gcv.hmin, vertices);
else
    v = NwSmooth(node, vin, radius, vertices);
end

end

function chanlocs = load_standard_channels()
% Load and configure standard 10-20 system channels
chanlocs = readlocs('Standard-10-20-Cap19.ced');

% Channel label standardization
new_labels = {'T3', 'T4', 'T5', 'T6'};
old_labels = {'T7', 'T8', 'P7', 'P8'};

for i = 1:numel(old_labels)
    idx = find(strcmpi({chanlocs.labels}, old_labels{i}));

    if ~isempty(idx)
        chanlocs(idx).labels = new_labels{i};
    end

end

end
