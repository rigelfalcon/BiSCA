% setup.m -- Add BiSCA toolbox paths to MATLAB search path.

addpath(genpath('./utility/'))
addpath(genpath('package/CNAR/'))
addpath(genpath('package/nufftkreg/'))
addpath(genpath('package/HoXiAlpha/'))
addpath(genpath('package/Pearson3/'))

% Vendored plotting deps (avoid requiring full installations)
if exist('./external/eeglab_min', 'dir') == 7
    addpath('./external/eeglab_min')
end
if exist('./external/Violinplot-Matlab', 'dir') == 7
    addpath('./external/Violinplot-Matlab')
end

set(groot, 'defaultFigureCloseRequestFcn', 'delete(gcf)');

% Allow running on newer NVIDIA GPUs with older CUDA libs.
try
    parallel.gpu.enableCUDAForwardCompatibility(true);
catch
end

% NOTE: do not call savepath here; it can hang in batch environments.
