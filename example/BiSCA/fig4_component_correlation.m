%% fig4_pearson_r
% Script:   fig4_component_correlation.m
% Figure:   Fig. 4c-f (component correlations)
% Purpose:  Calculate Pearson correlation coefficients for Fig. 4 panels.
%           Computes correlations between Full signal and Xi/Rho components
%           for both Gaussianity (median) and Linearity (max) statistics.
% Usage:    matlab -batch "run('./setup.m'); run('./example/BiSCA/fig4_component_correlation.m')"
%
% Output:   Prints Pearson r values for panels C-F to console and saves to file.
%

%% Setup
close all;
clearvars -except stat_bh;

% Resolve repo root
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

% Task name (same as fig4_ieeg_global_comp_bc.m)
taskname = 'Bxa15kFull2501051806';

%% Load data
if ~exist('stat_bh', 'var')
    load(fullfile(repo_root, 'data', 'result', 'stat_bh_summary.mat'))
end

%% Extract bicoherence magnitude (bc_T) for each component
% Gaussianity test uses 'median', Linearity test uses 'max'
%
% IMPORTANT: Fig.4C-F scatter plots use log10 scale on both axes.
% Therefore, Pearson r should be computed on log10-transformed data
% to match what the reader sees in the figure.

% Log10 transform (same as fig4_ieeg_global_comp_bc.m line 63)
op = @(x) log10(x);

% Full signal
bc_full_G = op([stat_bh.raw.median.bc_T]');  % Gaussianity (median-based)
bc_full_L = op([stat_bh.raw.max.bc_T]');     % Linearity (max-based)

% Xi component (aperiodic)
bc_xi_G = op([stat_bh.xi.median.bc_T]');
bc_xi_L = op([stat_bh.xi.max.bc_T]');

% Rho component (oscillatory) - stored as 'alpha' in data structure
bc_rho_G = op([stat_bh.alpha.median.bc_T]');
bc_rho_L = op([stat_bh.alpha.max.bc_T]');

%% Calculate Pearson correlations for Fig.4 panels C-F (iEEG)

% Panel C: Xi vs Full (Gaussianity)
[r_C, p_C] = corr(bc_full_G, bc_xi_G, 'Type', 'Pearson');

% Panel D: Rho vs Full (Gaussianity)
[r_D, p_D] = corr(bc_full_G, bc_rho_G, 'Type', 'Pearson');

% Panel E: Xi vs Full (Linearity)
[r_E, p_E] = corr(bc_full_L, bc_xi_L, 'Type', 'Pearson');

% Panel F: Rho vs Full (Linearity)
[r_F, p_F] = corr(bc_full_L, bc_rho_L, 'Type', 'Pearson');

%% Display results
fprintf('\n========================================\n');
fprintf('Pearson Correlation Coefficients\n');
fprintf('========================================\n');
fprintf('iEEG Data (N = %d channels)\n\n', length(bc_full_G));

fprintf('Panel C: Xi vs Full (Gaussianity, b^(G))\n');
fprintf('  Pearson r = %.4f, p = %.2e\n\n', r_C, p_C);

fprintf('Panel D: Rho vs Full (Gaussianity, b^(G))\n');
fprintf('  Pearson r = %.4f, p = %.2e\n\n', r_D, p_D);

fprintf('Panel E: Xi vs Full (Linearity, b^(L))\n');
fprintf('  Pearson r = %.4f, p = %.2e\n\n', r_E, p_E);

fprintf('Panel F: Rho vs Full (Linearity, b^(L))\n');
fprintf('  Pearson r = %.4f, p = %.2e\n\n', r_F, p_F);

%% Summary table
fprintf('Summary Table (for manuscript):\n');
fprintf('| Panel | Comparison | Test | Pearson r | p-value |\n');
fprintf('|-------|------------|------|-----------|----------|\n');
fprintf('| C | Xi vs Full | Gaussianity | %.3f | %.2e |\n', r_C, p_C);
fprintf('| D | Rho vs Full | Gaussianity | %.3f | %.2e |\n', r_D, p_D);
fprintf('| E | Xi vs Full | Linearity | %.3f | %.2e |\n', r_E, p_E);
fprintf('| F | Rho vs Full | Linearity | %.3f | %.2e |\n', r_F, p_F);
fprintf('\n');

%% Save results to file
outdir = fullfile(repo_root, 'result', 'pearson_r');
if ~exist(outdir, 'dir'), mkdir(outdir); end

results.taskname = taskname;
results.N = length(bc_full_G);
results.panel_C = struct('comparison', 'Xi_vs_Full', 'test', 'Gaussianity', 'r', r_C, 'p', p_C);
results.panel_D = struct('comparison', 'Rho_vs_Full', 'test', 'Gaussianity', 'r', r_D, 'p', p_D);
results.panel_E = struct('comparison', 'Xi_vs_Full', 'test', 'Linearity', 'r', r_E, 'p', p_E);
results.panel_F = struct('comparison', 'Rho_vs_Full', 'test', 'Linearity', 'r', r_F, 'p', p_F);

save(fullfile(outdir, 'pearson_results.mat'), 'results');
fprintf('Results saved to: %s\n', fullfile(outdir, 'pearson_results.mat'));
