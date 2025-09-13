function nar = nar_lpr(x, p, varargin)
    % NAR_LPR - Nonparametric Autoregression using Local Polynomial Regression
    % Input:
    %   x - Input time series (Tx1)
    % Output:
    %   nar - Model structure containing:
    %         X, Y - Hankel matrix inputs/outputs
    %         Y_smooth - Smoothed outputs
    %         hmin - Optimal bandwidth
    %         p - Autoregressive order
    %         H - Original Hankel matrix
    %         armodel - Best AR model
    if nargin < 2 || isempty(p)
        p_candidates = 1:200;
        armodel = fit_best_ar(x(:, 1), p_candidates, 'arburg');
        p = length(armodel.A);
    else

        armodel = [];
    end
    % Parse optional parameters
    pin = inputParser;
    addParameter(pin, 'order', 0, @isnumeric);
    addParameter(pin, 'isadaptive', false, @islogical)

    parse(pin, varargin{:});
    order = pin.Results.order;
    isadaptive = pin.Results.isadaptive;

    % Create Hankel matrix
    H = mvts2hankel(x, p);
    X = H(1:end - 1, :);
    Y = H(2:end, :);

    % Nadaraya-Watson smoothing
    if ~isadaptive
        % Bandwidth selection
        hlist = get_hList(100, [(max(H, [], 'all') - min(H, [], 'all')) / 100, ...
            (max(H, [], 'all') - min(H, [], 'all')) / 1], @logspace);
        [~, gcv_global] = NwSmoothGCVOriginalScaleLogGaussian(X, Y, hlist, X);
        hmin = gcv_global.hmin;

        [Y_smooth] = NwSmoothLogGaussianBatch(X, Y, hmin, X);
        kmin = [];
        h_i = [];
    else
        hmin = [];
        klist = unique(ceil(linspace(1, round(size(X, 1) / 20), 100))');
        [~, kmin, h_i, gcv_global] = NwSmoothMassCenterGCV(X, Y, klist, X);
        [Y_smooth] = NwSmoothMassCenterBatch(X, Y, kmin, X, h_i);

    end

    % Calculate residuals
    residuals = Y - Y_smooth;
    residual_cov = cov(residuals); % Full covariance matrix
    % Ensure positive definiteness for Cholesky decomposition
    try
        residual_chol = chol(residual_cov, 'lower');
    catch
        residual_cov = residual_cov + 1e-6 * eye(size(residual_cov));
        residual_chol = chol(residual_cov, 'lower');
    end

    % Store results in nar structure
    nar.X = X;
    nar.Y = Y;
    nar.order = order;
    nar.isadaptive = isadaptive;
    nar.hmin = hmin;
    nar.kmin = kmin;
    nar.h_i = h_i;
    nar.p = p;
    nar.armodel = armodel;
    nar.gcv = gcv_global;
    nar.residual_chol = residual_chol;
    nar.residual_cov = residual_cov;

end

function [yq, L, dbg] = NwSmoothLogGaussian(x, y, h, xq)
    % Log-Gaussian kernel smoother with numerical stability

    if nargin < 4 || isempty(xq)
        xq = x;
    end

    [n, dx] = size(x);
    [nq, ~] = size(xq);
    y = permute(y, [1, 3, 2]);
    % Dimension management
    if size(x, 2) > 1
        x = reshape(x, [n, 1, dx]);
    end
    if size(xq, 2) > 1
        xq = reshape(xq, [nq, 1, dx]);
    end

    h = permute(h, [dx + 2:-1:1]);
    D = x - permute(xq, [2, 1, 3]);

    % Log-Gaussian kernel computation
    log_kn = -0.5 * sum((D ./ h).^2, 3);

    % Numerical stability
    max_log = max(log_kn, [], 1);
    log_sum = log(sum(exp(log_kn - max_log), 1)) + max_log;
    weights = exp(log_kn - log_sum);

    % Prediction calculation
    yq = permute((sum(weights .* y, 1) ./ sum(weights, 1)), [2, 3, 1]);

    % Debug outputs
    dbg.s = sum(weights, 1);
    if nargout > 1
        L = weights;
    end
end

function [yq, dbg] = NwSmoothBatch(x, y, h, xq)
    % Process NW smoothing in batches to avoid memory overflow

    [n, dx] = size(x);
    nq = size(xq, 1);
    yq = zeros(nq, size(y, 2), 'like', y);
    x_exp = permute(x, [1, 3, 2]);
    xq_exp = permute(xq, [3, 1, 2]);
    y = permute(y, [1, 3, 2]);
    % Process query points in batches
    if isgpuarray(x)
        [bsize] = bytesize(1, 'b', single(1));
        [~, freemem] = test_gpu([], [], true);
    else
        [bsize] = bytesize(1, 'b', underlyingType(x(1)));
        [~, freemem] = test_memory();
    end
    num_free = max(1, floor(freemem / (bsize.b))); % heuristic

    num_per_batch = (num_free / (n * dx));
    num_batches = ceil(nq / num_per_batch);
    for i = 1:num_batches
        batch_idx = int64((i - 1) * num_per_batch + 1:min(i * num_per_batch, nq));
        xq_batch = xq_exp(:, batch_idx, :);

        D = x_exp - xq_batch;

        % Compute kernel weights and predictions
        Dkn = gaussian_kernel(D, h);
        sum_weights = sum(Dkn, 1);
        valid = sum_weights > eps('single');
        yq_batch = sum(Dkn .* y, 1) ./ (sum_weights + ~valid * eps('single'));
        yq(batch_idx, :) = gather(yq_batch); % Move to CPU if using GPU

        % Clear temporary GPU arrays
        clear D Dkn sum_weights yq_batch;
    end
    dbg.s = [];
end

function u = gaussian_kernel(u, b)
    u = u ./ b;
    u = (1 / sqrt(2 * pi)) * exp(-0.5 * sum((u .^ 2), 3));
end

function [yq, L, dbg] = NwSmoothLogGaussianBatch(x, y, h, xq)
    % Log-Gaussian kernel smoother using pdist2 and numerical stability
    % Inputs:
    %   x  - Input data [n × dx]
    %   y  - Response [n × dy]
    %   h  - Bandwidth [1 × dx] or scalar
    %   xq - Query points [nq × dx]
    % Outputs:
    %   yq - Predictions [nq × dy]
    %   L  - Weight matrix [nq × n]
    %   dbg - Debug info

    % Default query points
    if nargin < 4 || isempty(xq)
        xq = x;
    end

    % Dimension parameters
    [n, dx] = size(x);
    [nq, ~] = size(xq);
    dy = size(y, 2);

    % Bandwidth handling
    if isscalar(h)
        h = repmat(h, 1, dx); % Expand to [1 × dx]
    end

    % GPU memory check
    use_gpu = isgpuarray(x);
    if use_gpu
        [~, freemem] = test_gpu([], [], true);
        [elem_bytes] = bytesize(1, 'b', single(1));
    else
        [~, freemem] = test_memory();
        [elem_bytes] = bytesize(1, 'b', underlyingType(x(1)));
    end

    % Batch size calculation
    mem_per_batch = n * (dx + dy) * elem_bytes.b;
    batch_size = max(1, floor(0.8 * freemem / mem_per_batch));

    % Initialize outputs
    yq = zeros(nq, dy, 'like', y);
    if nargout > 1
        L = zeros(nq, n, 'like', y);
    end

    % Batch processing
    for bidx = 1:ceil(nq / batch_size)
        batch_range = (bidx - 1) * batch_size + 1:min(bidx * batch_size, nq);
        xq_batch = xq(batch_range, :);

        % Scaled distances [batch_size × n]
        scaled_x = x ./ h;
        scaled_q = xq_batch ./ h;
        D_sq = pdist2(scaled_q, scaled_x, 'squaredeuclidean');

        % Log-Gaussian kernel
        log_kn = -0.5 * D_sq;

        % Numerical stabilization
        max_log = max(log_kn, [], 2); % [batch_size × 1]
        log_weights = log_kn - max_log;
        log_sum = log(sum(exp(log_weights), 2)) + max_log;

        % Normalized weights
        weights = exp(log_kn - log_sum); % [batch_size × n]

        % Store outputs
        yq(batch_range, :) = weights * y;
        if nargout > 1
            L(batch_range, :) = weights;
        end
    end

    % Debug info
    if nargout > 2
        dbg.s = sum(weights, 2);
    end
end