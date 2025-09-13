%% fig4_component_metric_surfplots
% Script:   fig4_spatial_distribution.m
% Figure:   Fig. 4g-l (spatial distribution)
% Purpose:  Plot spatial distributions of component metrics used in Fig. 4.
% Usage:    matlab -batch "run('./setup.m'); run('./example/BiSCA/fig4_spatial_distribution.m')"
%

close all;
rng(0);

%% Configuration
% Resolve repo root (robust to current folder)
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

taskname = 'Bxa15kFull2501051806'; % Selected task
filepath = fullfile(repo_root, 'data', 'result', 'iEEG_bixialpha.mat');
figpath_fig3 = fullfile(repo_root, 'figure', 'fig3', filesep);
figpath_fig4 = fullfile(repo_root, 'figure', 'fig4', filesep);
test_folder(figpath_fig3);
test_folder(figpath_fig4);

opt.thfactor = 1; % Outlier removal threshold
opt.vfun = []; % No custom value function for scatter plots
opt.smooth_factor = 3;
%% Data Loading and Preparation
if ~exist('iEEG', 'var')
    load(filepath);
end

%%

if ~exist('bc', 'var')
    bixialpha = cat_struct_fields(iEEG.subjinfo, 'bixialpha', false, true);
    % Compute f_max_s and f_max_bc
    f_max_s = nan(length(bixialpha), 1);
    f_max_bc = nan(length(bixialpha), 1);
    f = bixialpha(1).f;
    frange_alpha = [3, 25];
    f_alpha = f(f >= frange_alpha(1) & f <= frange_alpha(2));

    for isubj = length(bixialpha):-1:1

        % Signal metrics
        [shat_xi.sum(isubj, :), shat_xi.mean(isubj, :), shat_xi.max(isubj, :)] = sum_mean_max(bixialpha(isubj).shat_xi);
        [shat_alpha.sum(isubj, :), shat_alpha.mean(isubj, :), shat_alpha.max(isubj, :)] = sum_mean_max(bixialpha(isubj).shat_alpha);
        shat_alpha_i = bixialpha(isubj).shat_alpha;
        shat_alpha_i = shat_alpha_i(f >= frange_alpha(1) & f <= frange_alpha(2), :);
        [~, i_max_s] = max(shat_alpha_i);
        f_max_s(isubj, :) = f_alpha(i_max_s);
        bc_alpha_i = bixialpha(isubj).bchat_alpha;
        bc_alpha_i = bc_alpha_i(f >= frange_alpha(1) & f <= frange_alpha(2), :);
        [~, i_max_bc] = max(abs(diag(bc_alpha_i)));
        f_max_bc(isubj, :) = f_alpha(i_max_bc);

        if isfield(bixialpha(isubj), 'para_alpha') && isfield(bixialpha(isubj).para_alpha, 'kernel')
            mus = bixialpha(isubj).para_alpha.kernel.para.mu;
        else
            mus = [];
        end
        bchat = bixialpha(isubj).bchat;
        F = griddedInterpolant({f, f}, bchat,"nearest");
        bc_mus_mus = F({mus, mus});
        bc_mus_mus = diag(bc_mus_mus);
        [~, i_max_bc_mus_mus] = max(abs(bc_mus_mus));
        f_max_bc_mus_mus(isubj, :) = mus(i_max_bc_mus_mus);

    end

end

%%

%%
if ~exist('stat_bh', 'var')
    filepath = fullfile(repo_root, 'data', 'result', 'stat_bh_summary.mat');
    load(filepath);
end

bc_G_xi = [stat_bh.xi.median.bc_T]';
bc_G_alpha = [stat_bh.alpha.median.bc_T]';
bc_L_xi = [stat_bh.xi.max.bc_T]';
bc_L_alpha = [stat_bh.alpha.max.bc_T]';

%% Main Plots: Scatter Surfaces and Histograms

% Mapping to manuscript panels:
% - Fig.3E: f_max_s (peak frequency of spectrum after removing Xi trend)
% - Fig.3F: f_max_bc (peak frequency of diagonal bicoherence after removing Xi trend)
% - Fig.4G: shat_xi.sum (Xi energy; sum of Xi component spectrum)
% - Fig.4H: shat_alpha.sum (Rho energy; sum of Rho component spectrum)
% - Fig.4I: bc_G_xi (Gaussianity statistic for Xi)
% - Fig.4J: bc_G_alpha (Gaussianity statistic for Rho)
% - Fig.4K: bc_L_xi (Nonlinearity statistic for Xi)
% - Fig.4L: bc_L_alpha (Nonlinearity statistic for Rho)
%
% Note: f_max_bc_mus_mus is an internal diagnostic (peak frequency on diagonal
% evaluated at kernel mu's); not shown in the manuscript figure.

% Log10 color mapping + unified color limits within same-quantity
% panel pairs.  G/H share energy clim; I/J share b^(G) clim; K/L share
% b^(L) clim.  Log10 compresses the dynamic range so that both xi and
% alpha panels show colour structure under a shared scale.
% Percentile [p5, p95] on the log10-transformed data avoids outlier stretch.
pctl = [5, 95];
energy_all = [shat_xi.sum; shat_alpha.sum];
energy_all = energy_all(~isnan(energy_all) & ~isinf(energy_all) & energy_all > 0);
bG_all     = [bc_G_xi; bc_G_alpha];
bG_all     = bG_all(~isnan(bG_all) & ~isinf(bG_all) & bG_all > 0);
bL_all     = [bc_L_xi; bc_L_alpha];
bL_all     = bL_all(~isnan(bL_all) & ~isinf(bL_all) & bL_all > 0);
clim_energy = prctile(log10(energy_all), pctl)';
clim_bG     = prctile(log10(bG_all), pctl)';
clim_bL     = prctile(log10(bL_all), pctl)';
fprintf('[fig4] clim_energy = [%.4f, %.4f] (log10, p%d-p%d)\n', clim_energy, pctl);
fprintf('[fig4] clim_bG     = [%.4f, %.4f] (log10, p%d-p%d)\n', clim_bG, pctl);
fprintf('[fig4] clim_bL     = [%.4f, %.4f] (log10, p%d-p%d)\n', clim_bL, pctl);

log10safe = @(x) log10(max(x, eps));  % guard against log10(0)

% Plot configs: {metric, quantities, op, climits, figdir, figname}
% Mapping to manuscript panels:
% - Fig.3E: f_max_s         -> fig3/fig3e_fmax_spectrum
% - Fig.3F: f_max_bc        -> fig3/fig3f_fmax_bicoherence
% - Fig.4G: shat_xi.sum     -> fig4/fig4g_shat_xi
% - Fig.4H: shat_alpha.sum  -> fig4/fig4h_shat_rho
% - Fig.4I: bc_G_xi         -> fig4/fig4i_bG_xi
% - Fig.4J: bc_G_alpha      -> fig4/fig4j_bG_rho
% - Fig.4K: bc_L_xi         -> fig4/fig4k_bL_xi
% - Fig.4L: bc_L_alpha      -> fig4/fig4l_bL_rho
plot_configs = {
    'f_max_s',     {''}, @(x)(x),  [],          figpath_fig3, 'fig3e_fmax_spectrum';
    'f_max_bc',    {''}, @(x)(x),  [],          figpath_fig3, 'fig3f_fmax_bicoherence';
    'shat_xi',     {'sum'}, log10safe, clim_energy, figpath_fig4, 'fig4g_shat_xi';
    'shat_alpha',  {'sum'}, log10safe, clim_energy, figpath_fig4, 'fig4h_shat_rho';
    'bc_G_xi',     {''}, log10safe, clim_bG,       figpath_fig4, 'fig4i_bG_xi';
    'bc_G_alpha',  {''}, log10safe, clim_bG,       figpath_fig4, 'fig4j_bG_rho';
    'bc_L_xi',     {''}, log10safe, clim_bL,       figpath_fig4, 'fig4k_bL_xi';
    'bc_L_alpha',  {''}, log10safe, clim_bL,       figpath_fig4, 'fig4l_bL_rho';
    };
close all;

for i = 1:size(plot_configs, 1)
    metric = plot_configs{i, 1};
    quantities = plot_configs{i, 2};
    op = plot_configs{i, 3};
    opt.climits = plot_configs{i, 4};
    figpath_i = plot_configs{i, 5};
    figname = plot_configs{i, 6};
    opt.vfun = [];

    for j = 1:length(quantities)
        quantity = quantities{j};

        if isempty(quantity)
            V = op(eval(metric));
        else
            V = op(abs(eval([metric '.' quantity])));
        end

        if all(isnan(V))
            V = zeros(size(V));
        end

        % Scatter surface plot
        figure;
        plot_scattersurf(iEEG, V, opt);

        save_fig(gcf, figpath_i, figname);

        % Histogram plot
        figure;
        histogram(rmoutliers(V, "quartiles", "ThresholdFactor", opt.thfactor), 100);
        ax = gca;
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;

        title(['Histogram of ', figname]);
    end

end

%%
function [x_sum, x_mean, x_max] = sum_mean_max(x)
x_sum = sum(x, "all", "omitnan");
x_mean = mean(x, "all", "omitnan");
x_max = max(x, [], "all", "omitnan");
end

function plot_scattersurf(iEEG, V, opt)
opt.reflective = false;
ChannelPos = cat_struct_fields(iEEG.subjinfo, 'ChannelPositionMidBst', false, true);

idx_bad = isnan(V) | isinf(V);
ChannelPos = ChannelPos(~idx_bad, :);
V = V(~idx_bad);

[~, idx_outlier] = rmoutliers(V, "quartiles", "ThresholdFactor", opt.thfactor);
ChannelPos = ChannelPos(~idx_outlier, :);
V = V(~idx_outlier);

clf;
opt.inflated = false;
opt.tess = iEEG.cortex.tess;
radius = max(max(iEEG.cortex.tess.Vertices, [], 1) - min(iEEG.cortex.tess.Vertices, [], 1)) * opt.smooth_factor;
surfaceplot5view([], [], V, ChannelPos, radius, true, false, opt);
end

function save_fig(hf, figpath, figname)
FontSize = 8;
Position = [1000, 375, 820, 512];

hf.Position = Position;
hA = gca;
hA.XAxis.FontSize = FontSize;
hA.YAxis.FontSize = FontSize;
dpi = 600;
colormap(plt.inferno);
exportgraphics(hf, [figpath, figname, '_', num2str(dpi), 'dpi', '.png'], 'ContentType', 'auto', 'Resolution', dpi);
if ispc; print(hf, [figpath, figname, '.emf'], '-dmeta'); end
end
