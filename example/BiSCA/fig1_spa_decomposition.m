%% fig1_diagram
% Script:   fig1_spa_decomposition.m
% Figure:   Fig. 1b (SPA decomposition)
% Purpose:  Run SPA-style fit and plotting for a single iEEG channel.
% Usage:    matlab -batch "run('./setup.m'); run('./example/BiSCA/fig1_spa_decomposition.m')"
%

close all; clearvars

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

force_export = false;
force_fit = false;
if exist('force_export', 'var') && ~isempty(force_export), force_export = true; end
if exist('force_fit', 'var') && ~isempty(force_fit), force_fit = true; end

figpath = fullfile(repo_root, 'figure', 'fig1', filesep);
outbase = fullfile(figpath, 'fig1f_spa_spectrum');
out_png = [outbase, '.png'];
out_svg = [outbase, '.svg'];


% Load demo channel data
demo_file = fullfile(repo_root, 'data', 'iEEG', 'iEEG_demo_channel.mat');
demo = load(demo_file);

ks = 0;
kh = 3;

% Cache the fitted BiXiAlpha object so we don't re-fit on every export.
cache_dir = fullfile(repo_root, 'data', 'result', 'figcache');
cachefile = fullfile(cache_dir, 'fig1f_spa_xialpha_subj32_ch3_ks0_kh3.mat');
test_folder(cachefile)

if exist(cachefile, 'file') == 2 && ~force_fit
    S = load(cachefile, 'xialpha');
    xialpha = S.xialpha;
    Fs = xialpha.Fs;
    fprintf('[fig1_SPA] loaded cached fit: %s\n', cachefile);
else
    x = double(demo.demo_signal);
    x = x(x ~= 0);
    Fs = demo.SamplingFrequency;

    hos = HOS(x, Fs);
    hos.window = 300;
    hos.nfft = 300;
    hos.frange = [0, 45];
    hos.overlap = 0.75;
    hos.normalization = 'haubrich'; % haubrich skewness
    hos.NW = 1.5; %1.5
    hos.method = 'pmtm'; % fft pmtm

    ks = 0;
    kh = 3; %3
    xialpha = BiXiAlpha(hos.f, hos.s, [], ks, kh, Fs);

    xialpha.method_alpha = 'findmax';
    xialpha.multiple_seg = false;
    xialpha.verbose = true;
    xialpha.para_alpha.no_zero_peak = true;
    xialpha.peak_relation = 'free';
    xialpha.regularization = 'elastic_constraint_mu_sigma_chi_xi';
    xialpha.lambda = [0; 0; 0; 0; 0]; %mu sigma startend noneg
    xialpha.peak_relation = 'free';
    xialpha.para_fit.StepTolerance = 1e-12;
    xialpha.para_fit.FunctionTolerance = 1e-12;
    xialpha.para_fit.OptimalityTolerance = 1e-12;
    xialpha.para_fit.FiniteDifferenceStepSize = 1e-8;
    xialpha.para_fit.MaxFunctionEvaluations = 20000;
    xialpha.para_fit.MaxIterations = 10000;

    xialpha.fit;

    save(cachefile, 'xialpha', 'ks', 'kh', '-v7.3');
    fprintf('[fig1_SPA] saved cached fit: %s\n', cachefile);
end


%% Fig.1f: SPA spectrum components
hos.f = xialpha.f;

figure(203)
clf
plot_s_shat_spa(xialpha)
drawnow

outbase = fullfile(figpath, 'fig1f_spa_spectrum');
test_folder([outbase, '.png'])
exportgraphics(gcf, [outbase, '.png'], 'ContentType', 'auto', 'Resolution', 600);
print(gcf, [outbase, '.svg'], '-dsvg');
if ispc
    print(gcf, [outbase, '.emf'], '-dmeta');
end
