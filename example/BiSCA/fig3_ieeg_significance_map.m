%% fig3_ieeg_sig
% Script:   fig3_ieeg_significance_map.m
% Figure:   Fig. 3c (iEEG significance map)
% Purpose:  Plot iEEG spatial distribution of Gaussianity/linearity test results.
% Usage:    matlab -batch "run('./setup.m'); run('./example/BiSCA/fig3_ieeg_significance_map.m')"
%

clearvars -except iEEG stat_bh

% close all
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

if ~exist('iEEG', 'var')
    iEEG = load(fullfile(repo_root, 'data', 'sEEG', 'sEEG_metadata.mat'));
end

iEEG = remove_useless_ieeg_fields(iEEG);

if ~exist('stat_bh', 'var')
    load(fullfile(repo_root, 'data', 'result', 'stat_bh_summary.mat'))
end

%%
p_threshold = 0.001;

comp = 'raw'; %raw alpha alpha
type = ['maxmedian'];
quantity = 'h'; %zscore


%%

figure(101)
clf
Nch_subj = cat_struct_fields(iEEG.subjinfo, 'NumChannel', false, true);
idx_remove = Nch_subj == 1;
iEEG.subjinfo(idx_remove) = [];
Nch_subj(idx_remove) = [];
iEEG.Nsubj = length(Nch_subj);
idx_subj_ch = gen_idx_subj_ch(Nch_subj);
test_method = 'stat_bh'; %stat stat_btstp_ar stat_btstp_kernel_ar stat_btstp_comp stat_btstp_kernel_comp

v = zeros(length(idx_subj_ch), 1);
varname = [comp, '_', type, '_', quantity];

for idx = length(idx_subj_ch):-1:1
    isubj = idx_subj_ch(idx, 1);
    ich = idx_subj_ch(idx, 2);

    if ich == Nch_subj(isubj)
        iEEG.subjinfo(isubj).(varname) = [];
    end

    iEEG.subjinfo(isubj).(varname)(ich, :) = 1+bin2dec([num2str(eval(test_method).(comp).('max')(idx).(quantity)),num2str(eval(test_method).(comp).('median')(idx).(quantity))]);
    v(idx) = iEEG.subjinfo(isubj).(varname)(ich, :);
end

opt.isbinary = false;
opt.reflective=true; %true false
opt.cmp=colormap(flip(jet(4)));
ieeg_scatterplot5view(iEEG, varname, opt);
figpath = fullfile(repo_root, 'figure', 'fig3', filesep);
figname = fullfile(figpath, 'fig3_iEEG_test_colormap_');
test_folder(figpath);
dpi=900;
exportgraphics(gcf, [figname, num2str(dpi), 'dpi', '.png'], 'ContentType', 'auto', 'Resolution', dpi);
%%
figure(111)
plot_class_pie(v,{'L+G','L+NG','NL+G','NL+NG'});
figname = fullfile(figpath, 'fig3_iEEG_test_pie_');
exportgraphics(gcf, [figname, '.eps'], 'ContentType', 'vector');
