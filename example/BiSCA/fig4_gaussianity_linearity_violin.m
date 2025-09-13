%% fig4_ieeg_global_comp_bc
% Script:   fig4_gaussianity_linearity_violin.m
% Figure:   Fig. 4a-b (Gaussianity/linearity violin plots)
% Purpose:  Plot component-wise Gaussianity/nonlinearity statistics used in Fig. 4.
% Usage:    matlab -batch "run('./setup.m'); run('./example/BiSCA/fig4_gaussianity_linearity_violin.m')"
%

%% Common Setup
taskname = 'Bxa15kFull2501051806';
close all;

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

% Load data
if ~exist('stat_bh', 'var')
    load(fullfile(repo_root, 'data', 'result', 'stat_bh_summary.mat'))
end
%%
% Visualization parameters
xycmp = plt.viridis(3);

% Create color matrix with group-specific hues
color_matrix = [xycmp(1, :) * 0.4; xycmp(1, :);
                xycmp(2, :) * 0.6; xycmp(2, :);
                xycmp(3, :) * 0.6; xycmp(3, :)];

labels = {'${b}^{(L)}_{full}$', ...
              '${b}^{(G)}_{full}$'; ...
              '${b}^{(L)}_{\xi}$', ...
              '${b}^{(G)}_{\xi}$'; ...
              '${b}^{(L)}_{\rho}$', ...
          '${b}^{(G)}_{\rho}$'};

FontSize = 16;
linewidth = 4;
dpi = 900;
p_threshold = 0.001;
Position = [200, 200, 500, 500];

%% Figure 501: Component-wise comparisons (supports Fig. 4C-F)
% This section generates 2x2 comparison grids between the full signal and its
% components, consistent with the Fig.4 scatterplot narrative (C-F). Here the
% plotted quantity is bc_T (bicoherence magnitude) on log10 scale.
figure(501)
clf;
f = gcf;
f.Position = [277 199.6667 566 521];

% Data organization
test_method = 'stat_bh';
quantity = 'bc_T'; % Changed from zscore to T

types = {'max', 'median'};
op = @(x) log10(x); %log10

for t = 1:2
    type = types{t};
    val.raw.(type).(quantity) = op([eval(test_method).raw.(type).(quantity)]');
    val.xi.(type).(quantity) = op([eval(test_method).xi.(type).(quantity)]');
    val.alpha.(type).(quantity) = op([eval(test_method).alpha.(type).(quantity)]');
end

% Threshold calculation (T-value specific)
Ns = eval(test_method).raw.max.Ns;

bc_th_max = mean([eval(test_method).raw.(types{1}).bc_th]);
bc_th_median = mean([eval(test_method).raw.(types{2}).bc_th]);
threshold = op([bc_th_max, bc_th_median]);

% Create plots
tiledlayout(2, 2, 'TileSpacing', 'loose');
plot_comparison_grid(val, threshold, quantity, color_matrix, 14, linewidth);

% Save figure (intermediate grid plot used for Fig.4C-F exploration)
% Version 1: without colorbar
figname = fullfile(repo_root, 'figure', 'fig4', 'fig4cf_bct_fdr');
fprintf('[fig4] Saving C-F to: %s\n', figname);
test_folder(fullfile(repo_root, 'figure', 'fig4', filesep));
exportgraphics(f, [figname, '_', num2str(dpi), 'dpi.png'], 'Resolution', dpi);
if ispc; print(f, [figname, '.emf'], '-dmeta'); end
fprintf('[fig4] C-F saved.\n');

%%

% Create a figure for the legend
figure(401); clf
hold on;

% Plot small patches or lines for each color
for i = 1:size(color_matrix, 1)
    % Plot a small line segment or patch for the color
    plot([1 1.5], [6 - i 6 - i], 'LineWidth', 6, 'Color', color_matrix(i, :));
    % Add the corresponding label with LaTeX interpreter
    text(1.7, 6 - i, labels{i}, 'Interpreter', 'latex', 'FontSize', FontSize);
end

% Customize the figure
axis([0 5 0 7]); % Adjust axis limits to fit the legend
set(gca, 'XTick', [], 'YTick', []); % Remove axes ticks
box off; % Remove the box around the plot
set(gca, 'Visible', 'off'); % Hide axes entirely
hold off;

%% Figure 502: Violin plots (supports Fig. 4A-B)
% Left subplot: Gaussianity statistics (b^(G)) for full / Xi / Rho (Fig.4A)
% Right subplot: Nonlinearity statistics (b^(L)) for full / Xi / Rho (Fig.4B)
% Both use log10 scale and show threshold lines.

test_method = 'stat_bh';
figure(502)
clf;
f = gcf;
f.Position = [459 776.3333 520 400];
FontSize = 14;

% Data preparation with logarithmic transformation
epsilon = 1e-5;
violin_data = cell(3, 2); % [groups x metrics]

% Max values
violin_data{1, 1} = [eval(test_method).raw.max.bc_T]' + epsilon;
% Mean values
violin_data{1, 2} = [eval(test_method).raw.median.bc_T]' + epsilon;

% Xi values
violin_data{2, 1} = [eval(test_method).xi.max.bc_T]' + epsilon;
violin_data{2, 2} = [eval(test_method).xi.median.bc_T]' + epsilon;

% Alpha values
violin_data{3, 1} = [eval(test_method).alpha.max.bc_T]' + epsilon;
violin_data{3, 2} = [eval(test_method).alpha.median.bc_T]' + epsilon;

% Apply logarithmic transformation to data
log_data = cellfun(@(x) log10(x), violin_data, 'UniformOutput', false);

% Create plot positions and labels
labels = {'${b}^{(G)}$', ...
              '${b}^{(G)}_{\xi}$', ...
              '${b}^{(G)}_{\rho}$', ...
              '${b}^{(L)}$', ...
              '${b}^{(L)}_{\xi}$', ...
              '${b}^{(L)}_{\rho}$'};

% Separate L and G data
log_data = log_data';
log_data = log_data(:);
combined_log_data = cat(2, log_data{:});

% Split data into G and L tests
g_data = combined_log_data(:, 2:2:end); % G tests
l_data = combined_log_data(:, 1:2:end); % L tests

% Create subplots
subplot(1, 2, 1)
violinplot(g_data, labels(1:3), ...
    'Width', 0.6, ...
    'Orientation', 'vertical', ...
    'ViolinColor', color_matrix(2:2:end, :), ...
    'ShowMean', true, ...
    'MarkerSize', 2, ...
    'MedianMarkerSize', 40, ...
    'BoxWidth', 0.03, ...
    'EdgeColor', [0.5 0.5 0.5]);

subplot(1, 2, 2)
violinplot(l_data, labels(4:6), ...
    'Width', 0.6, ...
    'Orientation', 'vertical', ...
    'ViolinColor', color_matrix(1:2:end, :), ...
    'ShowMean', true, ...
    'MarkerSize', 2, ...
    'MedianMarkerSize', 40, ...
    'BoxWidth', 0.03, ...
    'EdgeColor', [0.5 0.5 0.5]);

% Configure both subplots
for i = 1:2
    subplot(1, 2, i)
    min_log = floor(min([g_data(:); l_data(:)]));
    max_log = ceil(max([g_data(:); l_data(:)]));
    yticks = min_log:max_log;
    yticklabels = arrayfun(@(x) sprintf('$10^{%d}$', x), yticks, 'UniformOutput', false);

    % Add threshold lines
    hold on;
    bc_th_max = mean([eval(test_method).raw.(types{1}).bc_th]);
    bc_th_median = mean([eval(test_method).raw.(types{2}).bc_th]);
    thresholds = log10([bc_th_max, bc_th_median]);

    % Plot only the relevant threshold for each subplot
    if i == 1 % Gaussian Tests
        th_y = thresholds(2);
        yline(th_y, '--k', 'LineWidth', 2);
    else % Linear Tests
        th_y = thresholds(1);
        yline(th_y, '--k', 'LineWidth', 2);
    end

    % Annotate above/below threshold at top/bottom edges
    yl = [min_log, max_log];
    xl = xlim;
    if i == 1
        lbl_above = 'Non-Gaussian';
        lbl_below = 'Consistent with Gaussian';
    else
        lbl_above = 'Non-linear';
        lbl_below = 'Consistent with linear';
    end
    text(xl(2), yl(2), lbl_above, ...
        'FontSize', FontSize - 5, 'Color', [0.3 0.3 0.3], ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
    text(xl(2), yl(1), lbl_below, ...
        'FontSize', FontSize - 5, 'Color', [0.3 0.3 0.3], ...
        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

    % Configure axes
    if i == 1
        xlabels = labels(1:3);
    else
        xlabels = labels(4:6);
    end

    set(gca, 'XTick', 1:3, ...
        'XTickLabel', xlabels, ...
        'YTick', yticks, ...
        'YTickLabel', yticklabels, ...
        'YLim', [min_log max_log], ...
        'FontSize', FontSize, ...
        'TickLabelInterpreter', 'latex');
    ylabel('${b}$ ($log_{10}$ scale)', 'Interpreter', 'latex');
    grid on;
end

% Save figure (violin plot for Fig.4A-B)
figname = fullfile(repo_root, 'figure', 'fig4', 'fig4ab_violin');
fprintf('[fig4] Saving violin to: %s\n', figname);
exportgraphics(f, [figname, '.png'], 'Resolution', 600);
if ispc; print(f, [figname, '.emf'], '-dmeta'); end
fprintf('[fig4] Violin saved.\n');

%% Helper Functions

function export_figure(f, figname, taskname, dpi)
% Standardized figure export
test_folder(figname);
exportgraphics(f, [figname, taskname, '_', num2str(dpi), 'dpi.png'], 'Resolution', dpi);
if ispc; print(f, [figname, taskname, '.emf'], '-dmeta'); end
end

function x = dlog(x)
% Modified log transform
x = sign(x) .* log10(1 + abs(x));
end

function h = hist_log_marg(x, y, xlims, ylims, xycmp)

if nargin < 3 || isempty(xlims)
    xlims = [];
end

if nargin < 3 || isempty(ylims)
    ylims = [];
end

h = plot_histogram_with_marginals(x, y, xlims, ylims, xycmp);
end

function plot_comparison_grid(val, thresholds, quantity, xycmp, FontSize, linewidth)
% PLOT_COMPARISON_GRID Creates 2x2 grid of T-value comparisons
%
% Disable per-tile density colorbar; use a single shared one at bottom.
setappdata(gcf, 'show_density_colorbar', false);

% Create tiled layout
t = tiledlayout(2, 2);
t.TileSpacing = 'loose';

lims = {[min(val.raw.max.(quantity)), max((val.raw.max.(quantity))); ...
             min(val.xi.max.(quantity)), max(val.xi.max.(quantity))]; ...

    [min(val.raw.median.(quantity)), max((val.raw.median.(quantity))); ...
     min(val.xi.median.(quantity)), max(val.xi.median.(quantity))]; ...

    [min(val.raw.max.(quantity)), max((val.raw.max.(quantity))); ...
     min(val.alpha.max.(quantity)), max(val.alpha.max.(quantity))]; ...

    [min(val.raw.median.(quantity)), max((val.raw.median.(quantity))); ...
     min(val.alpha.median.(quantity)), max(val.alpha.median.(quantity))]};

% Plot pairs (Fig.4C-F): full vs components
panel_labels = {'C', 'E', 'D', 'F'};  % Order matches tile layout (row-major)

plot_config = {
    % Fig.4C
    {'median', 'xi', thresholds(2), 'G', 'G', '\xi', xycmp([2, 4], :), lims{2}}
    % Fig.4E
    {'max', 'xi', thresholds(1), 'L', 'L', '\xi', xycmp([1, 3], :), lims{1}}
    % Fig.4D
    {'median', 'alpha', thresholds(2), 'G', 'G', '\rho', xycmp([2, 6], :), lims{4}}
    % Fig.4F
    {'max', 'alpha', thresholds(1), 'L', 'L', '\rho', xycmp([1, 5], :), lims{3}}

    };

for i = 1:4
    nexttile;
    cfg = plot_config{i};
    type = cfg{1};
    comp = cfg{2};
    thresh = cfg{3};
    x_sup = cfg{4};
    y_sup = cfg{5};
    comp_name = cfg{6};
    color = cfg{7}; % Get color from config
    lim = cfg{8};
    xlims = lim(1, :);
    ylims = lim(2, :);

    xlims(1) = min(xlims(1), thresh) - 0.1;
    ylims(1) = min(ylims(1), thresh) - 0.1;
    xlims(2) = max(xlims(2), thresh) + 0.1;
    ylims(2) = max(ylims(2), thresh) + 0.1;

    % Get data pair
    x_data = val.raw.(type).(quantity);
    y_data = val.(comp).(type).(quantity);

    % Plot histogram with marginals
    hist_log_marg(x_data, y_data, xlims, ylims, color);

    % Add threshold lines
    xline(thresh, 'LineWidth', linewidth, 'Color', [0.5 0.5 0.5]);
    yline(thresh, 'LineWidth', linewidth, 'Color', [0.5 0.5 0.5]);

    % Set labels
    xlabel(['${b}^{(', x_sup, ')}$ ($log_{10}$ scale)'], 'Interpreter', 'latex', 'FontSize', FontSize);
    ylabel(['${b}^{(', y_sup, ')}_{', comp_name, '}$ ($log_{10}$ scale)'], 'Interpreter', 'latex', 'FontSize', FontSize);

    % Configure axes
    ax = gca;
    ax.XAxis.FontSize = FontSize;
    ax.YAxis.FontSize = FontSize;

    grid on;
    box on;
end

% Unify density colorbar across ALL 4 panels (C-F).
tiles = t.Children;
ax_tiles = flipud(tiles(arrayfun(@(h) isa(h, 'matlab.graphics.axis.Axes'), tiles)));
cl_all = [min(arrayfun(@(a) a.CLim(1), ax_tiles)), max(arrayfun(@(a) a.CLim(2), ax_tiles))];
for k = 1:numel(ax_tiles)
    clim(ax_tiles(k), cl_all);
end
fprintf('[fig4_CF] unified density clim ALL = [%.2f, %.2f]\n', cl_all);

end
