function [yq, nw_gcv] = NwSmoothGCVOriginalScaleLogGaussian(x, yr, hList, xq, N, israndomize)
    % Unified Nadaraya-Watson smoother with GCV using log-Gaussian kernel

    use_gpu = true;
    [Tx, dx] = size(x);
    dh = size(hList, 1);
    [Ty, dy] = size(yr);

    % Default parameter handling
    if nargin < 5 || isempty(N), N = ceil(Tx^(1 / dx)); end
    if nargin < 4 || isempty(xq), xq = get_ndgrid_scatter(x, 'list', N); end
    if nargin < 6, israndomize = []; end

    % Memory management
    if use_gpu
        [x, yr, hList] = convert_to_gpuarray(x, yr, hList);
        % Avoid test_gpu() since it may call gpuDevice(count) and re-select
        % devices, which can invalidate existing gpuArrays on some setups.
        memCheckFcn = @() local_gpu_available_mem();
    else
        memCheckFcn = @() test_memory();
    end
    if size(hList, 2) > 1 && size(hList, 2) ~= size(x, 2)
        error(['column of hList must match to the column of x']);
    end
    freemem = memCheckFcn();
    Txr = set_subsample_size(Tx, dx, freemem, use_gpu, israndomize);

    % Subsampling if needed
    if Txr < Tx
        idx = sort(randsample(Tx, Txr, false));
        x = x(idx, :);
        yr = yr(idx, :);
    end

    % Main processing loop
    [gcvm, pdof, df, mse] = process_bandwidths_log(x, yr, hList, dx, dy, dh, Txr);

    % GCV results packaging
    [nw_gcv] = package_results(x, gcvm, pdof, df, mse, hList, dx);

    [yq] = NwSmoothLogGaussianBatch(x, yr, nw_gcv.hmin, xq);

    disp(['[nw_smooth_log]: [', num2str(dx), 'd] [', num2str(size(x, 1)), ' samples] regression']);
end

function [gcvm, pdof, df, mse] = process_bandwidths_log(x, yr, hList, dx, dy, dh, Txr)
    yhat = zeros(Txr, dy, dh, 'like', yr);
    [pdof, gcvm, df, mse] = deal(zeros(dh, 1));

    for i = 1:dh
        h = hList(i, :);
        [yq, diagL] = NwSmoothLogGaussianBatch(x, yr, h, x);
        yhat(:, :, i) = yq;
        % Diagnostics calculation
        pdof(i) = 1 - mean(diagL);
        df(i) = sum(diagL);
        mse(i) = mean((yr - yq).^2, 'all');
        gcvm(i) = mse(i) / pdof(i)^2;
    end
end

function Txr = set_subsample_size(Tx, dx, freemem, use_gpu, israndomize)
    % Determine optimal subsample size based on memory constraints
    % Inputs:
    %   Tx        - Original sample size [scalar]
    %   dx        - Data dimensions [scalar]
    %   freemem   - Available memory in bytes [scalar]
    %   use_gpu   - GPU flag [logical]
    %   israndomize - Manual override flag [logical | empty]
    % Output:
    %   Txr       - Reduced sample size [scalar]

    % Determine element byte size
    if use_gpu
        bytes_per_element = 4; % Single-precision for GPU
    else
        bytes_per_element = 8; % Double-precision for CPU
    end

    % Original code's q-factor (1 = conservative, 1.5 = aggressive)
    q_factor = 1;

    if isempty(israndomize)
        % Auto-detect based on memory requirements
        required_memory = Tx^2 * dx * bytes_per_element;
        israndomize = required_memory > (q_factor * freemem);
    end

    if israndomize
        % Calculate max affordable sample size
        Txr = floor(sqrt((q_factor * freemem) / (dx * bytes_per_element)));
        Txr = min(Txr, Tx); % Never exceed original size
    else
        Txr = Tx; % Use full dataset
    end
end

function freemem = local_gpu_available_mem()
    d = gpuDevice;
    freemem = d.AvailableMemory;
end

function [nw_gcv] = package_results(x, gcvm, pdof, df, mse, hList, dx)
    [gcvm_min, idmin] = min(gcvm);

    nw_gcv = struct(...
        'idmin', idmin, ...
        'gcvm', gcvm, ...
        'pdof', pdof, ...
        'df', df, ...
        'mse', mse, ...
        'hmin', hList(idmin, :), ...
        'pdof_min', pdof(idmin));

    disp(['[nw_smooth]: [', num2str(dx), 'd] [', num2str(size(x, 1)), ' samples] regression'])
    disp([' - selected bandwidth - [', num2str(find(ismember(hList, nw_gcv.hmin.', 'rows'))), 'th] h'])
    disp([' - normalized scale h =', num2str(nw_gcv.hmin)])
    disp([' - original scale h =', num2str(nw_gcv.hmin .* std(x))])
    disp([' - gcvm =', num2str(gcvm_min)])
    disp([' - mse =', num2str(mse(idmin))])
end