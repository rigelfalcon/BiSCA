%% fig3_ieeg_freq_bct
% Script:   fig3_ieeg_bifrequency_test.m
% Figure:   Fig. 3d (iEEG bifrequency test)
% Purpose:  Plot frequency summaries related to bicoherence test statistics.
% Usage:    matlab -batch "run('./setup.m'); run('./example/BiSCA/fig3_ieeg_bifrequency_test.m')"
%

close all
clearvars -except stat_bh

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
stat_all=stat_bh;
alpha=stat_all.raw.median(1).p_threshold;

% Extract channel and bixialpha data
Nch_subj = cat_struct_fields(iEEG.subjinfo, 'NumChannel', false, true);
idx_subj_ch = gen_idx_subj_ch(Nch_subj);
bixialpha_all = cat_struct_fields(iEEG.subjinfo, 'bixialpha', false, true);
%%
stat_all.max_abs_bc_f=[];
stat_all.max_abs_bc=[];
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



%% Plotting with FDR Correction and Symlog Scale on Y-Axis (One-Tailed Test) using Robust LWP
% Create the plot
figure(111); clf; hold on;
set(gcf, "Position", [790.3333, 554.3333, 1024, 512]);
colormap(plt.viridis)
cmap = colormap(plt.viridis(256));

xlims=[0,25];
ylims=[-1.5,0.5];


type_stat = 'bc';
type_metric = 'median';

% Process median data
x = reshape(stat_all.max_abs_bc_f', [], 1);
y = [stat_all.raw.(type_metric).bc_T]';
y = reshape((y .* ones(1, size(stat_all.max_abs_bc_f,2)))', [], 1);
[~,idx_median]=rmoutliers(y,"quartiles", "ThresholdFactor",3);
x = x(~idx_median);
y = y(~idx_median);

% Mask for valid data
mask = ~isnan(y) & ~isinf(y); % No need to exclude zero for symlog
x_median = x(mask);
y_median = y(mask);

% Process max data
type_metric = 'max';
x = reshape(stat_all.max_abs_bc_f', [], 1);
y = [stat_all.raw.(type_metric).bc_T]';
y = reshape((y .* ones(1, size(stat_all.max_abs_bc_f,2)))', [], 1);
[~,idx_max]=rmoutliers(y,"quartiles", "ThresholdFactor",3);
x = x(~idx_max);
y = y(~idx_max);

% Mask for valid data
mask = ~isnan(y) & ~isinf(y);
x_max = x(mask);
y_max = y(mask);

% Combine x data for histogram
x_all = [x_median; x_max];

% FDR Correction for Z-Scores (One-tailed test)
% Extract p-values corresponding to z-scores
p_vals = [stat_all.raw.median.p]'; % Assuming p-values are stored here
num_tests = length(p_vals); % 1772 channels

% Apply FDR (Benjamini-Hochberg procedure)
[sorted_p, sort_idx] = sort(p_vals);
fdr_threshold = (1:num_tests)' * alpha / num_tests; % FDR level = alpha
fdr_reject = sorted_p <= fdr_threshold;
last_reject = find(fdr_reject, 1, 'last');
if ~isempty(last_reject)
    p_fdr = sorted_p(last_reject); % Critical p-value
else
    p_fdr = 0; % No significant results
end

bc_th_median = [stat_all.raw.median.bc_th]; % One-tailed threshold for positive z-scores
bc_th_median = bc_th_median(~idx_median);
bc_th_median=mean(log10(bc_th_median),'all','omitmissing');

bc_th_max = [stat_all.raw.max.bc_th]; % One-tailed threshold for positive z-scores
bc_th_max = bc_th_max(~idx_max);
bc_th_max=mean(log10(bc_th_max),'all','omitmissing');

% Apply symlog transformation to y-data
y_median_scale = log10(y_median);
y_max_scale = log10(y_max);

% Main plot
main_ax = gca;
xlim(xlims)
ylim(ylims);

scatter(x_median+0.5, y_median_scale, 10, 'filled', 'MarkerFaceColor', cmap(1, :), 'MarkerFaceAlpha', 0.1, ...
    'MarkerEdgeColor', cmap(1, :),'MarkerEdgeAlpha',0.2, 'DisplayName', '$b_G$ of samples');

scatter(x_max, y_max_scale, 10, 'filled', 'MarkerFaceColor', cmap(end, :), 'MarkerFaceAlpha', 0.1, ...
    'MarkerEdgeColor', cmap(end, :),'MarkerEdgeAlpha',0.2, 'DisplayName', '$b_L$ of samples');

FontSize = 12;
hA = gca;


hA.XAxis.FontSize = FontSize;
hA.YAxis.FontSize = FontSize;

% Add trend lines using robust LWP
x_vals = linspace(min(x_median(x_median > 0)), max(x_median), 1000)';
hlist=get_hList(100,[2,10],@logspace);
y_vals_median=NwSmoothLogGaussianBatch(x_median,y_median_scale,2,x_vals);
plot(x_vals, y_vals_median, 'LineWidth', 2, 'Color', max(min(cmap(1, :)*1.5,[1,1,1]),[0,0,0]), 'DisplayName', 'Trend of $b_G$');

x_vals_max = linspace(min(x_max(x_max > 0)), max(x_max), 1000)';
y_vals_max=NwSmoothLogGaussianBatch(x_max,y_max_scale,2,x_vals_max);
plot(x_vals_max, y_vals_max, 'LineWidth', 2, 'Color', max(min(cmap(end, :)*0.8,[1,1,1]),[0,0,0]), 'DisplayName', 'Trend of $b_L$');

% Format main plot
xlabel('Frequency (Hz)');
ylabel('$b$', 'Interpreter', 'latex'); % Updated label for symlog scale
grid on;
hLegend = legend('Location', 'northeastoutside', 'Interpreter', 'latex', 'FontSize', FontSize);
uistack(hLegend, 'top');

% Median threshold line: solid blue line
yline(bc_th_median, 'LineWidth', 3, ...
     'Color', [0, 0, 1], ...
     'LineStyle', '--', ...
     'DisplayName', ['FDR p=' num2str(alpha) ' for $b_G$']);

% Max threshold line: dashed red line
yline(bc_th_max, 'LineWidth', 3, ...
     'Color', [1, 0, 0], ...
     'LineStyle', '--', ...
     'DisplayName', ['FDR p=' num2str(alpha) ' for $b_L$']);

% Adjust main axis position
main_ax_pos = main_ax.Position;
main_ax.Position = [main_ax_pos(1), main_ax_pos(2), main_ax_pos(3) * 0.85, main_ax_pos(4) * 0.85];

% Top histogram (x density)
top_ax = axes('Position', [main_ax_pos(1), main_ax_pos(2) + main_ax_pos(4) * 0.85, main_ax_pos(3) * 0.85, main_ax_pos(4) * 0.15]);
histogram(top_ax, x_median, 'FaceColor', cmap(128, :), 'EdgeColor', 'none', 'BinWidth', 2,'Normalization','pdf');
top_ax.XLim = main_ax.XLim;
top_ax.YAxis.Visible = 'off';
top_ax.XAxis.Visible = 'off';
uistack(top_ax, 'bottom');

% Right histogram (y density)
right_ax = axes('Position', [main_ax_pos(1) + main_ax_pos(3) * 0.85, main_ax_pos(2), main_ax_pos(3) * 0.15, main_ax_pos(4) * 0.85]);
hold(right_ax, 'on');
histogram(right_ax, y_median_scale, 'Orientation', 'horizontal', 'FaceColor', cmap(1, :), 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'BinWidth', 0.1,'Normalization','pdf');
histogram(right_ax, y_max_scale, 'Orientation', 'horizontal', 'FaceColor', cmap(end, :), 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'BinWidth', 0.1,'Normalization','pdf');
ylim(ylims);
right_ax.XAxis.Visible = 'off';
right_ax.YAxis.Visible = 'off';
uistack(right_ax, 'bottom');

% Save figure
figpath = fullfile(repo_root, 'figure', 'fig3', filesep);
figname = ['fig3_iEEG_pfp_vs_bct_', type_stat];
save_fig(gcf, figpath, figname);

function save_fig(hf, figpath, figname)

f = gcf;
dpi = 600;
test_folder(figpath);
exportgraphics(hf, [figpath, figname, '_', num2str(dpi), 'dpi', '.png'], 'ContentType', 'auto', 'Resolution', dpi);
end
