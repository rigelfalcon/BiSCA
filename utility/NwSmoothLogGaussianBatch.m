function [yq, diagL, L, dbg] = NwSmoothLogGaussianBatch(x, y, h, xq)
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
    % Check if xq is x for diagonal computation
    computeDiag = nargout > 1 && isequaln(xq, x);
    if computeDiag
        diagL = zeros(nq, 1, 'like', y);
    end
    % GPU memory check
    use_gpu = isgpuarray(x);
    if use_gpu
        % Avoid test_gpu() since it may call gpuDevice(count) and re-select
        % devices, which can invalidate existing gpuArrays on some setups.
        d = gpuDevice;
        freemem = d.AvailableMemory;
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
        currentSize = numel(batch_range);

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
        if computeDiag
            idx = sub2ind([currentSize, n], 1:currentSize, batch_range);
            diagL(batch_range, :) = weights(idx);
        else
            diagL = [];
        end

        if nargout > 2
            L(batch_range, :) = weights;
        end
    end

    % Debug info
    if nargout > 2
        dbg.s = sum(weights, 2);
    end
end