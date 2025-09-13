function [yq, h_i, diagL, L, dbg] = NwSmoothMassCenterBatch(x, y, k, xq, h_i)
    % Mass-center kernel smoother with adaptive bandwidths (batched for memory efficiency)
    % Inputs:
    %   x  - Data points (n x dx matrix)
    %   y  - Response values (n x dy matrix)
    %   k  - Number of neighbors for bandwidth estimation (default: sqrt(n))
    %   xq - Query points (nq x dx matrix, default: x)
    % Outputs:
    %   yq - Smoothed values at query points (nq x dy)
    %   L  - Weight matrix (nq x n) if requested
    %   dbg- Debug information with bandwidths and sums

    % Handle default inputs
    if nargin < 4 || isempty(xq)
        xq = x;
    end

    if nargin < 3 || isempty(k)
        k = ceil(sqrt(size(x, 1)));
    end

    % Dimension parameters
    [n, dx] = size(x);
    [nq, ~] = size(xq);
    dy = size(y, 2);

    % Centralized GPU Data Management
    % This function is the single source of truth for device state.
    % All incoming data is assumed to be on the CPU, and this block
    % will move it to the GPU if required.

    % Determine if GPU should be used based on the first input array.
    % This is a proxy for the 'use_gpu' option passed from the top level.
    use_gpu = isgpuarray(x);

    if use_gpu
        % Ensure all necessary variables are on the GPU.
        [x, y, k, xq, h_i] = convert_to_gpuarray(x, y, k, xq, h_i);
    end

    if nargin < 5 || isempty(h_i)
        % This part should ideally not be called when using GCV, but is kept for standalone use.
        % Perform knnsearch on CPU for stability, then move h_i to GPU if needed.
        [~, dist] = knnsearch(gather(x), gather(x), 'K', k + 1);
        h_i_cpu = sqrt((dist(:, end) .^ 2) ./ dx);
        if use_gpu
            h_i = gpuArray(h_i_cpu);
        else
            h_i = h_i_cpu;
        end
    end

    % Precompute normalization factor for Gaussian kernel
    norm_factor = 1;

    % Determine memory constraints and batch size
    if use_gpu
        % test_gpu() can call gpuDevice(count) and re-select device, which may
        % invalidate existing gpuArrays during forward-compat compilation.
        % Query current device memory without re-selecting.
        d = gpuDevice;
        freemem = d.AvailableMemory;
        [elem_bytes] = bytesize(1, 'b', single(1));
    else
        [~, freemem] = test_memory();
        [elem_bytes] = bytesize(1, 'b', underlyingType(x(1)));
    end

    % Memory required per batch (n * batch_size * dx elements for 3D array)
    mem_per_batch = n * dx * elem_bytes.b;
    batch_size = max(1, floor(0.8 * freemem / mem_per_batch));

    % Initialize outputs
    yq = zeros(nq, dy, 'like', y);

    if nargout > 2
        diagL = zeros(nq, 1, 'like', y); % only work when nq == n
    end

    if nargout > 3
        L = zeros(nq, n, 'like', y);
    end

    if nargout > 4
        dbg.s = zeros(nq, 1, 'like', y);
        dbg.h_i = h_i;
    end

    % Process in batches
    for bidx = 1:ceil(nq / batch_size)
        batch_range = (bidx - 1) * batch_size + 1:min(bidx * batch_size, nq);
        xq_batch = xq(batch_range, :);
        current_batch_size = numel(batch_range);

        % Compute scaled differences and squared distances
        D_sq = pdist2(x, xq_batch, 'squaredeuclidean'); %squaredeuclidean fastsquaredeuclidean for cpu
        sum_sq = D_sq ./ (h_i .^ 2);

        % Compute Gaussian weights
        weights = norm_factor .* exp(-0.5 * sum_sq);
        normalized_weights = weights ./ sum(weights, 1);
        % Store results
        yq(batch_range, :) = normalized_weights' * y;

        if nargout > 2

            if nq == n
                idx = sub2ind(size(normalized_weights), batch_range, 1:current_batch_size);
                diagL(batch_range, :) = normalized_weights(idx);
            else
                diagL = [];
            end

        end

        if nargout > 3
            L(batch_range, :) = normalized_weights';
        end

        if nargout > 4
            dbg.s(batch_range) = sum(weights, 1);
        end

    end

end