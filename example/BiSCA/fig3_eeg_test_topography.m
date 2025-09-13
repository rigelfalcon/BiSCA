%% fig3_eeg_spatial_freq_statistic
% Script:   fig3_eeg_test_topography.m
% Figure:   Fig. 3a-b (EEG test topographies)
% Purpose:  Generate spatial and frequency summary statistics for rsEEG.
% Usage:    matlab -batch "run('./setup.m'); run('./example/BiSCA/fig3_eeg_test_topography.m')"
%

close all
clearvars -except stat_all eeg_summary

% Get repo root (match iEEG pattern)
script_path = mfilename('fullpath');
[script_dir, ~, ~] = fileparts(script_path);
repo_root = fileparts(fileparts(script_dir));

% Guard EEGLAB dependency: expect setup.m to add it.
if ~exist('topoplot', 'file')
    warning('topoplot not found. Ensure EEGLAB is on path via setup.m');
    warning('Skipping figure generation because EEG topographies require topoplot.');
    return;
end

%%
close all;

% Load pre-aggregated EEG summary
summary_file = fullfile(repo_root, 'data', 'result', 'HarMNqEEG', 'eeg_stat_summary.mat');
if ~exist(summary_file, 'file')
    warning('EEG summary not found: %s', summary_file);
    warning('Run tools/strip_data_for_release.m to generate, then re-run.');
    return;
end

if ~exist('eeg_summary', 'var')
    eeg_summary = load(summary_file);
end
stat_all = eeg_summary.stat_all;

% Reconstruct age info in stat_all for compatibility
if isfield(eeg_summary, 'age')
    nch = numel(stat_all.raw.max);
    age_vec = eeg_summary.age;
    if numel(age_vec) == nch
        for ki = 1:nch
            stat_all.raw.max(ki).info.age = age_vec(ki);
        end
    end
end

%% proportion
HG = [stat_all.raw.median.h];
HL = [stat_all.raw.max.h];
Total = length(HG);
G = sum(double(~HG));
NG = sum(double(HG));
L = sum(double(~HL));
NL = sum(double(HL));

p_G = G / Total;
p_L = L / Total;
p_NG = NG / Total;
p_NL = NL / Total;
disp(['p_G: ', num2str(p_G), ' p_L: ', num2str(p_L), ' p_NG: ', num2str(p_NG), ' p_NL: ', num2str(p_NL)]);

%%
GL = ~HG & ~HL;
NGL = HG & ~HL;
GNL = ~HG & HL;
NGNL = HG & HL;

disp(['GL: ', num2str(sum(GL)), ' NGL: ', num2str(sum(NGL)), ' GNL: ', num2str(sum(GNL)), ' NGNL: ', num2str(sum(NGNL))]);

p_GL = sum(GL) / Total;
p_NGL = sum(NGL) / Total;
p_GNL = sum(GNL) / Total;
p_NGNL = sum(NGNL) / Total;
disp(['p_GL: ', num2str(p_GL), ' p_NGL: ', num2str(p_NGL), ' p_GNL: ', num2str(p_GNL), ' p_NGNL: ', num2str(p_NGNL)]);

%%

% subplot a: topography and probability of HG and HL
chan_names = {'Fp1' 'Fp2' 'F3' 'F4' 'C3' 'C4' 'P3' 'P4' 'O1' ...
              'O2' 'F7' 'F8' 'T3' 'T4' 'T5' 'T6' 'Fz' 'Cz'}';

type_stat = 'P';
type_metric = 'max';

y = GL';
y = y([stat_all.raw.(type_metric).Ns_orig] > 140, :);
y = reshape((y)', 18, [])';
p_GL = y;
p_GL = (sum(p_GL, 1) ./ size(p_GL, 1))';

% Process max data
y = NGL';
y = y([stat_all.raw.(type_metric).Ns_orig] > 140, :);
y = reshape((y)', 18, [])';
p_NGL = y;
p_NGL = (sum(p_NGL, 1) ./ size(p_NGL, 1))';

y = GNL';
y = y([stat_all.raw.(type_metric).Ns_orig] > 140, :);
y = reshape((y)', 18, [])';
p_GNL = y;
p_GNL = (sum(p_GNL, 1) ./ size(p_GNL, 1))';

y = NGNL';
y = y([stat_all.raw.(type_metric).Ns_orig] > 140, :);
y = reshape((y)', 18, [])';
p_NGNL = y;
p_NGNL = (sum(p_NGNL, 1) ./ size(p_NGNL, 1))';

%%

op = @(x) x;
% Unified color scale across 4 topo panels (EEG internal)
clims = [min([p_GL; p_NGL; p_GNL; p_NGNL], [], 'all'), max([p_GL; p_NGL; p_GNL; p_NGNL], [], 'all')];
fprintf('[fig3_eeg] topo clims = [%.4f, %.4f]\n', clims(1), clims(2));
figure(301); clf
POS=[458.3333333333333, 476.3333333333333, 200, 180];

figure(301); clf
set(gcf, "Position", POS);

y = p_GL;
plot_topography(y, 'ChannelLabel', chan_names, 'SigThreshold', [], 'CBTitle', '$\mathit{p}_{GL}$')
clim(clims)
colorbar

filename = fullfile(repo_root, 'figure', 'fig3', 'fig3_rsEEG_topo_1.png');
test_folder(filename);
exportgraphics(gcf, filename, 'Resolution', 600);
if ispc; print(gcf, strrep(filename, '.png', '.emf'), '-dmeta'); end

figure(302); clf
set(gcf, "Position", POS);
y = p_NGL;
plot_topography(y, 'ChannelLabel', chan_names, 'SigThreshold', [], 'CBTitle', '$\mathit{p}_{NGL}$')
clim(clims)
colorbar

filename = fullfile(repo_root, 'figure', 'fig3', 'fig3_rsEEG_topo_2.png');
test_folder(filename);
exportgraphics(gcf, filename, 'Resolution', 600);
if ispc; print(gcf, strrep(filename, '.png', '.emf'), '-dmeta'); end

figure(303); clf
set(gcf, "Position", POS);
y = p_GNL;
plot_topography(y, 'ChannelLabel', chan_names, 'SigThreshold', [], 'CBTitle', '$\mathit{p}_{GNL}$')
clim(clims)
colorbar
filename = fullfile(repo_root, 'figure', 'fig3', 'fig3_rsEEG_topo_3.png');
test_folder(filename);
exportgraphics(gcf, filename, 'Resolution', 600);
if ispc; print(gcf, strrep(filename, '.png', '.emf'), '-dmeta'); end


figure(304); clf
set(gcf, "Position", POS);
y = p_NGNL;
plot_topography(y, 'ChannelLabel', chan_names, 'SigThreshold', [], 'CBTitle', '$\mathit{p}_{NGNL}$')
clim(clims)
colorbar
filename = fullfile(repo_root, 'figure', 'fig3', 'fig3_rsEEG_topo_4.png');
test_folder(filename);
exportgraphics(gcf, filename, 'Resolution', 600);
if ispc; print(gcf, strrep(filename, '.png', '.emf'), '-dmeta'); end

%%
% Load pre-computed frequency probability maps
pf = eeg_summary.p_freq;
p_freq_GL   = pf.p_freq_GL;
p_freq_NGL  = pf.p_freq_NGL;
p_freq_GNL  = pf.p_freq_GNL;
p_freq_NGNL = pf.p_freq_NGNL;
p_freq_NG   = pf.p_freq_NG;
p_freq_NL   = pf.p_freq_NL;
N = size(p_freq_GL, 1);
f = linspace(eeg_summary.frange(1), eeg_summary.frange(end), N);
[fx, fy] = meshgrid(f, f);

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

%%
figure(401)
set(gcf, 'Position', [197.6667, 991, 450, 350]);
% Unified linear color limits across all 4 bifrequency panels (EEG internal)
linear_data_max = max([p_freq_GL(:); p_freq_NGL(:); p_freq_GNL(:); p_freq_NGNL(:)], [], 'omitmissing');
fprintf('[fig3_eeg] linear data max = %.4f\n', linear_data_max);
clims = [0, linear_data_max];

tiledlayout(2, 2)
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

%save figure
filename = fullfile(repo_root, 'figure', 'fig3', 'fig3_rsEEG_freq_stat.png');
test_folder(filename);
exportgraphics(gcf, filename, 'Resolution', 600);
if ispc; print(gcf, strrep(filename, '.png', '.emf'), '-dmeta'); end

%%

%%
function plot_topography(channel_data, varargin)
% PLOT_TOPGRAPHY Create EEG channel topography plot with log color scaling
%   plot_topography(channel_data, 'Param', value, ...)
%
%   Inputs:
%   channel_data - Vector of values per channel (natural scale values)
%
%   Optional Parameters:
%   'Colormap'    - Custom colormap (default: inferno)
%   'CLim'        - Color axis limits in natural scale (default: [min max] of data)
%   'LogScale'    - Use log color scale (default: true)
%   'SigThreshold'- Significance threshold to mark on colorbar (natural scale)
%   'SavePath'    - Directory to save figure
%   'FileName'    - Base filename for saving
%   'Visible'     - Figure visibility ('on' or 'off')
%   'ChannelLabel'- Cell array of channel labels to include
%   'CBTitle'     - Colorbar title (default: 'Value')

%% Parameter parsing
p = inputParser;
addRequired(p, 'channel_data', @isnumeric);
addParameter(p, 'Colormap', (plt.inferno(256)), @(x) isnumeric(x) && size(x, 2) == 3);
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
% Requires EEGLAB topoplot on MATLAB path.
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

%% Helper function for channel configuration
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
