%% Nonlinear Autoregressive Local Polynomial Regression (NAR-LPR) for iEEG Data Analysis
% This script demonstrates nonlinear and non-Gaussian dynamics analysis using NAR-LPR
% on intracranial EEG (iEEG) data, including phase space reconstruction and potential landscape visualization.

%% Initialize environment
close all
clearvars -except iEEG % Keep iEEG data if already loaded
rng(0) % Set random seed for reproducibility

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

% Define initial parameters
N = 11600; % Initial length of time series

%% Load iEEG data if not already in workspace
if ~exist('iEEG', 'var')
    filepath = fullfile(repo_root, 'data', 'iEEG', 'iEEG1772FemTposHmInvScMleLog.mat');
    iEEG = load(filepath);
    iEEG = remove_useless_ieeg_fields(iEEG); % Clean up unnecessary fields
end

% Select specific subject and channel for analysis
[isubj, ich] = deal(94, 8); % Subject 94, channel 8

% Extract data segment and prepare for processing
x = double(iEEG.subjinfo(isubj).Data_W(1:N, ich));
x = x(501:2500, :); % Select specific time window
x = gpuArray(single(x)); % Convert to GPU array for faster computation

% Plot raw signal and update N to match actual data length
figure(1);
clf;
plot(x);
title('Raw iEEG Signal');
xlabel('Sample');
ylabel('Amplitude');
N = size(x, 1);

%% Set up figure export parameters
figpath = fullfile(repo_root, 'figure', 'nar_lpr_sym/'); % Directory for saving figures
dpi = 600; % Resolution for exported figures


%% Create time vector and plot time series
% Convert sample indices to time in seconds
t = (0:length(x) - 1)' / iEEG.SamplingFrequency;

% Plot time series with proper time axis
figure(2);
clf;
plot(t, x);
xlabel('Time (s)');
ylabel('Amplitude');

%% Fit NAR-LPR model
% Create nonlinear autoregressive local polynomial regression model
% Using order 23 with adaptive fitting
nar = nar_lpr(x, 23, 'isadaptive', true);

%% Initial NAR-LPR forecasting
% Find the state with maximum energy as initial condition
[~, idx_max] = max(sum(abs(nar.X).^2, 2));
init_stat = nar.X(idx_max, :);

% Generate forecast and visualize in 2D phase space
figure(101);
clf;
[y_pred, low_dim_traj_2d] = nar_lpr_forecast(nar, init_stat, N, ...
    'vis', true, ...
    'vis_train', true, ...
    'vis_dim', 2, ...
    'nburning', 100);

%% Prune NAR-LPR model to focus on equilibrium dynamics
% Screen equilibrium keys with threshold 1e-2
nar_pruned = nar_lpr_screen_equilibrium_keys(nar, y_pred, 'threshold', 1e-2);
nar_pruned.issym = true;
%% Generate forecast with pruned model
% Find maximum energy state in pruned model
[~, idx_max] = max(sum(abs(nar_pruned.X).^2, 2));
init_stat = nar_pruned.X(idx_max, :);

% Generate deterministic forecast (sigma=0) with pruned model
figure(101);
clf;
[y_pred_pruned, low_dim_traj_2d_pruned] = nar_lpr_forecast(nar_pruned, init_stat, N, ...
    'vis', true, ...
    'vis_train', true, ...
    'vis_dim', 2, ...
    'nburning', 0, ...
    'sigma', 0);
%% Visualize potential landscape in 2D
% Estimate potential landscape using PCA projection
sigma = 1; % Smoothing parameter for density estimation
num_steps = 200; % Number of grid points for potential estimation
result = estimate_potential_pca(nar_pruned, num_steps, 100, sigma, 1000);
%%
% Plot the potential landscape with trajectory overlay
figure(102);
clf;
plot_potential_pca(result, low_dim_traj_2d_pruned(1:num_steps * 2, :));
set(gca, 'SortMethod', 'depth'); % Proper depth sorting for 3D visualization
view(2); % Set to 2D top-down view

% Set figure size for better visualization
set(gcf, "Position", [485 300 820 700]);

% Save the figure
figname = [mfilename, '2d_potential_map'];
test_folder([figpath, figname, '_', num2str(dpi), 'dpi', '.png']);
exportgraphics(gcf, [figpath, figname, '_', num2str(dpi), 'dpi', '.png'], ...
    'ContentType', 'auto', 'Resolution', dpi);

%% Visualize 3D phase space portrait - Raw data
figure(201);
clf;
% Generate 3D phase space visualization with noise (sigma=1)
nar_lpr_forecast(nar_pruned, init_stat, N, ...
    'vis', true, ...
    'vis_train', true, ...
    'vis_forecast', false, ...
    'vis_dim', 3, ...
    'sigma', 1);

% Set optimal viewing angle
view(-215.4869, 37.8437);

% Save the figure
figname = [mfilename, '3d_phase_space_raw'];
test_folder([figpath, figname, '_', num2str(dpi), 'dpi', '.png']);
exportgraphics(gcf, [figpath, figname, '_', num2str(dpi), 'dpi', '.png'], ...
    'ContentType', 'auto', 'Resolution', dpi);

%% Visualize 3D phase space portrait - With noise
figure(201);
clf;
% Generate 3D phase space visualization with noise (sigma=1)
nar_lpr_forecast(nar_pruned, init_stat, N, ...
    'vis', true, ...
    'vis_train', true, ...
    'vis_dim', 3, ...
    'sigma', 1);

% Set optimal viewing angle
view(-215.4869, 37.8437);

% Save the figure
figname = [mfilename, '3d_phase_space_sim_noise'];
test_folder([figpath, figname, '_', num2str(dpi), 'dpi', '.png']);
exportgraphics(gcf, [figpath, figname, '_', num2str(dpi), 'dpi', '.png'], ...
    'ContentType', 'auto', 'Resolution', dpi);

%% Simulate with custom noise innovation
figure(202);
clf;

% Set burning period for initialization
nburning = 100;

% Generate custom Gaussian noise innovation scaled by model residual variance
innovation = ar_noise(11600 + nburning, nar_pruned.p, 'gaussian');

% Generate forecast with custom noise innovation
[y_pred_noise, low_dim_traj_3d_noise] = nar_lpr_forecast(nar_pruned, init_stat, 11600, ...
    'vis', true, ...
    'vis_train', true, ...
    'vis_dim', 3, ...
    'sigma', 1, ...
    'innovation_type', 'input', ...
    'innovation', innovation, ...
    'nburning', nburning);

%% Simulate free-running deterministic system (no noise)
figure(203);
clf;

% Generate forecast with no noise (sigma=0)
[y_pred_free, low_dim_traj_3d_free] = nar_lpr_forecast(nar_pruned, init_stat, 11600, ...
    'vis', true, ...
    'vis_train', true, ...
    'vis_dim', 3, ...
    'sigma', 0, ...
    'nburning', nburning);

% Set optimal viewing angle
view(-215.4869, 37.8437);

% Save the figure
figname = [mfilename, '3d_phase_space_sim_free'];
test_folder([figpath, figname, '_', num2str(dpi), 'dpi', '.png']);
exportgraphics(gcf, [figpath, figname, '_', num2str(dpi), 'dpi', '.png'], ...
    'ContentType', 'auto', 'Resolution', dpi);


%% Compare original and simulated time series
figure(501);
clf;

% Set figure size for better visualization
set(gcf, "Position", [860 450 640 240]);

% Extract portion of data for comparison
x_plot = x(1500:N);
t = (0:length(x_plot) - 1)' / iEEG.SamplingFrequency;

% Create timetable for stacked plot comparison
dataTable = timetable(seconds(t), ...
    gather(x_plot), ...
    y_pred_free(1500:N, 1), ...
    'VariableNames', {'Original', 'FreeModel'});

% Create stacked plot to compare time series
stackedplot(dataTable, 'FontSize', 12, 'Color', [0, 0, 0]);

% Save the figure
figname = [mfilename, 'timeseries'];
test_folder([figpath, figname, '_', num2str(dpi), 'dpi', '.png']);
exportgraphics(gcf, [figpath, figname, '_', num2str(dpi), 'dpi', '.png'], ...
    'ContentType', 'auto', 'Resolution', dpi);


%% Higher-Order Spectra (HOS) analysis - Raw data
% Initialize HOS object with raw data and sampling frequency
hos = HOS(x(:, 1), 200); % 200 Hz sampling frequency

% Configure HOS parameters
hos.window = 300; % Window size for spectral estimation
hos.nfft = 300; % FFT size
hos.frange = [0, 50]; % Frequency range of interest (Hz)
hos.overlap = 0.75; % Window overlap ratio
hos.normalization = 'haubrich'; % Normalization method for skewness
hos.NW = 2.5; % Time-bandwidth product for multitaper
hos.method = 'pmtm'; % Use multitaper method for spectral estimation
hos.tf_return_res = false; % Don't return residuals
hos.tf_return_seg = false; % Don't return segments
hos.log_var_s = true; % Use logarithmic scale for variance

% Compute bicoherence
hos.bc;

% Plot power spectrum
figure(401);
clf;
plot_s(hos.f, hos.s);

% Save power spectrum figure
figname = [mfilename, 'spectrum_raw'];
test_folder([figpath, figname, '_', num2str(dpi), 'dpi', '.png']);
exportgraphics(gcf, [figpath, figname, '_', num2str(dpi), 'dpi', '.png'], ...
    'ContentType', 'auto', 'Resolution', dpi);

% Plot bicoherence
figure(411);
clf;
plot_bc(hos.fx, hos.fy, hos.bc, @abs);
view(2);

% Set figure size and save
set(gcf, 'Position', [400 200 300 300]);
figname = [mfilename, 'bicoherence_raw'];
test_folder([figpath, figname, '_', num2str(dpi), 'dpi', '.png']);
exportgraphics(gcf, [figpath, figname, '_', num2str(dpi), 'dpi', '.png'], ...
    'ContentType', 'auto', 'Resolution', dpi);

figure(431);
clf;
Ns_eff = get_hos_Ns_eff(x(:, 1), hos, 100, false, true); %
visualize_bc_stat_max_median_symlog_log(hos.bc, Ns_eff, 0.001)

%% Higher-Order Spectra (HOS) analysis - Noise model
% Initialize HOS object with noise model data
hos = HOS(y_pred_noise(:, 1), 200); % 200 Hz sampling frequency

% Configure HOS parameters (same as for raw data)
hos.window = 300; % Window size for spectral estimation
hos.nfft = 300; % FFT size
hos.frange = [0, 50]; % Frequency range of interest (Hz)
hos.overlap = 0.75; % Window overlap ratio
hos.normalization = 'haubrich'; % Normalization method for skewness
hos.NW = 2.5; % Time-bandwidth product for multitaper
hos.method = 'pmtm'; % Use multitaper method for spectral estimation
hos.tf_return_res = false; % Don't return residuals
hos.tf_return_seg = false; % Don't return segments
hos.log_var_s = true; % Use logarithmic scale for variance

% Compute bicoherence
hos.bc;

% Plot power spectrum
figure(402);
clf;
plot_s(hos.f, hos.s);

% Save power spectrum figure
figname = [mfilename, 'spectrum_noise'];
test_folder([figpath, figname, '_', num2str(dpi), 'dpi', '.png']);
exportgraphics(gcf, [figpath, figname, '_', num2str(dpi), 'dpi', '.png'], ...
    'ContentType', 'auto', 'Resolution', dpi);

% Plot bicoherence
figure(412);
clf;
plot_bc(hos.fx, hos.fy, hos.bc, @abs);
view(2);

% Set figure size and save
set(gcf, 'Position', [400 200 300 300]);
figname = [mfilename, 'bicoherence_noise'];
test_folder([figpath, figname, '_', num2str(dpi), 'dpi', '.png']);
exportgraphics(gcf, [figpath, figname, '_', num2str(dpi), 'dpi', '.png'], ...
    'ContentType', 'auto', 'Resolution', dpi);

figure(432);
clf;
Ns_eff = get_hos_Ns_eff(y_pred_noise(:, 1), hos, 100, false, true); %
visualize_bc_stat_max_median_symlog_log(hos.bc, Ns_eff, 0.001)


%% Higher-Order Spectra (HOS) analysis - Free-running model
% Initialize HOS object with free-running model data
hos = HOS(y_pred_free(:, 1), 200); % 200 Hz sampling frequency

% Configure HOS parameters (same as for previous analyses)
hos.window = 300; % Window size for spectral estimation
hos.nfft = 300; % FFT size
hos.frange = [0, 50]; % Frequency range of interest (Hz)
hos.overlap = 0.75; % Window overlap ratio
hos.normalization = 'haubrich'; % Normalization method for skewness
hos.NW = 2.5; % Time-bandwidth product for multitaper
hos.method = 'pmtm'; % Use multitaper method for spectral estimation
hos.tf_return_res = false; % Don't return residuals
hos.tf_return_seg = false; % Don't return segments
hos.log_var_s = true; % Use logarithmic scale for variance

% Compute bicoherence
hos.bc;

% Plot power spectrum
figure(403);
clf;
plot_s(hos.f, hos.s);

% Save power spectrum figure
figname = [mfilename, 'spectrum_free'];
test_folder([figpath, figname, '_', num2str(dpi), 'dpi', '.png']);
exportgraphics(gcf, [figpath, figname, '_', num2str(dpi), 'dpi', '.png'], ...
    'ContentType', 'auto', 'Resolution', dpi);

% Plot bicoherence
figure(413);
clf;
plot_bc(hos.fx, hos.fy, hos.bc, @abs);
view(2);

% Set figure size and save
set(gcf, 'Position', [400 200 300 300]);
figname = [mfilename, 'bicoherence_free'];
test_folder([figpath, figname, '_', num2str(dpi), 'dpi', '.png']);
exportgraphics(gcf, [figpath, figname, '_', num2str(dpi), 'dpi', '.png'], ...
    'ContentType', 'auto', 'Resolution', dpi);

figure(433);
clf;
Ns_eff = get_hos_Ns_eff(y_pred_free(:, 1), hos, 100, false, true); %
visualize_bc_stat_max_median_symlog_log(hos.bc, Ns_eff, 0.001)


%% Helper functions

% Function to plot power spectrum with consistent formatting
function plot_s(f, s)
    % Set figure size
    set(gcf, 'Position', [400 200 300 200]);

    % Plot power spectrum in dB scale
    plot(f, 10 * log10(s), '-', 'LineWidth', 1, 'Color', [0, 0, 0]);
    hold on;

    % Configure plot appearance
    grid off;
    xlim([0, 50]);
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');

    % Set font sizes
    FontSize = 12;
    hA = gca;
    hA.XAxis.FontSize = FontSize;
    hA.YAxis.FontSize = FontSize;

    % Set colormap and tick marks
    colormap(hA, plt.viridis);
    dtick = 5;
    freqTicks = 0:dtick:max(f(:));
    xticks(hA, freqTicks);
    box off;
    set(gca, 'Layer', 'top');
    hold off;
end