%% fig1_diagram
% Script:   fig1_bisca_model_fit.m
% Figure:   Fig. 1c-d (spectrum + bicoherence model fit)
% Purpose:  Generate BiSCA-related Fig.1 panels (excluding SPA panel).
% Usage:    matlab -batch "run('./setup.m'); run('./example/BiSCA/fig1_bisca_model_fit.m')"
%

close all

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

% ---- Fast path: if exports exist, skip everything ----
force_export = false;
if exist('force_export', 'var') && ~isempty(force_export)
    force_export = true;
end

figpath = fullfile(repo_root, 'figure', 'fig1', filesep);
out_e_svg = fullfile(figpath, 'fig1e_spectrum_components.svg');
out_e_png = fullfile(figpath, 'fig1e_spectrum_components_600dpi.png');
out_d_svg = fullfile(figpath, 'fig1d_bicoherence_fit.svg');
out_d_png = fullfile(figpath, 'fig1d_bicoherence_fit_600dpi.png');

if exist(out_e_svg, 'file') == 2 && exist(out_e_png, 'file') == 2 && ...
        exist(out_d_svg, 'file') == 2 && exist(out_d_png, 'file') == 2 && ~force_export
    fprintf('[fig1_BiSCA] outputs exist, skip export: %s\n', figpath);
    return
end

if ~exist('iEEG', 'var')
    iEEG = load(fullfile(repo_root, 'data', 'sEEG', 'sEEG_metadata.mat'));
end

iEEG = remove_useless_ieeg_fields(iEEG);

% ---- Config (override in caller workspace before running) ----
if ~exist('taskname', 'var') || isempty(taskname)
    taskname = 'Bxa15kFull2501051806';
end

if ~exist('idx_ch', 'var') || isempty(idx_ch)
    idx_ch = 579; % demo global channel index
end

Nch_subj = cat_struct_fields(iEEG.subjinfo, 'NumChannel', false, true);
idx_subj_ch = gen_idx_subj_ch(Nch_subj);
isubjich = num2cell(idx_subj_ch(idx_ch, :));
[isubj, ich] = deal(isubjich{:});

if ~exist('fig1_result_file', 'var')
    fig1_result_file = '';
end

%%
savepath = fullfile(repo_root, 'data', 'result', 'iEEG', 'HoXiAlpha', ...
    ['W', filesep, taskname, '_cpct'], ...
    ['iEEG_', num2str(isubj), '_', num2str(ich), '.mat']);

if exist(savepath, 'file') == 2
    load(savepath, 'bixialpha');
else
    error('Missing required Fig.1 result file: %s', savepath);
end

iEEG.subjinfo(isubj).bixialpha(ich, :) = bixialpha;
test_folder(figpath)
figpath_aux = fullfile(figpath, 'aux', filesep);
test_folder(figpath_aux);


%% Figure 1e: BiSCA spectrum components
hos.f = bixialpha.f;
figure(203)
clf
plot_s_shat_struct(bixialpha)

ax = gca;
fprintf('[fig1_BiSCA] figure(203) children=%d\n', numel(ax.Children));
lgd = legend(ax);
if isgraphics(lgd)
    try
        fprintf('[fig1_BiSCA] figure(203) legend items=%d\n', numel(lgd.String));
    catch
    end
end

drawnow
outbase = fullfile(figpath, 'fig1e_spectrum_components');
print(gcf, [outbase, '.svg'], '-dsvg');

dpi = 600;
exportgraphics(gcf, [outbase, '_', num2str(dpi), 'dpi', '.png'], ...
    'ContentType', 'auto', 'Resolution', dpi);

if ispc
    print(gcf, [outbase, '.emf'], '-dmeta');
end


%% Figure 1d: BiSCA bicoherence (with harmonic guide lines)
figure(303)
clf
plot_bc_bchat_fit_scatter_struct(bixialpha)

outbase = fullfile(figpath, 'fig1d_bicoherence_fit');
dpi = 600;

exportgraphics(gcf, [outbase, '_', num2str(dpi), 'dpi', '.png'], 'ContentType', 'auto', 'Resolution', dpi);
print(gcf, [outbase, '.svg'], '-dsvg');
if ispc
    print(gcf, [outbase, '.emf'], '-dmeta');
end
