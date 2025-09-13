function [f] = kde_batch(x, xq, h)
    % KDE with pdist2 and batch processing for memory efficiency
    % Inputs:
    %   x  - Data points [n x d]
    %   xq - Query points [nq x d]
    %   h  - Bandwidth [scalar or 1 x d]
    % Output:
    %   f  - Density estimates [nq x 1]
    [x, xq, h] = convert_to_gpuarray(x, xq, h);

    [n, d] = size(x);
    nq = size(xq, 1);

    % Memory configuration
    use_gpu = true;

    if use_gpu
        % Avoid test_gpu() since it may call gpuDevice(count) and re-select
        % devices, which can invalidate existing gpuArrays on some setups.
        d0 = gpuDevice;
        free_mem = d0.AvailableMemory;
        [elem_bytes] = bytesize(1, 'b', single(1));
    else
        [~, free_mem] = test_memory();
        [elem_bytes] = bytesize(1, 'b', underlyingType(x(1)));
    end

    % Batch size calculation (conservative estimate)
    mem_per_point = n * elem_bytes.b; % Distance vector per query point
    batch_size = max(1, floor(0.8 * free_mem / mem_per_point));

    % Bandwidth handling
    if isscalar(h)
        h = repmat(h, 1, d);
    end

    % Precompute normalization constant
    log_norm_const = -d / 2 * log(2 * pi) - sum(log(h));

    % Batch processing
    f = zeros(nq, 1, 'like', x);
    for bidx = 1:ceil(nq / batch_size)
        % Get batch indices
        start_idx = (bidx - 1) * batch_size + 1;
        end_idx = min(bidx * batch_size, nq);
        batch_q = xq(start_idx:end_idx, :);

        % Compute squared distances [nq_batch x n]
        D_sq = pdist2(batch_q ./ h, x ./ h, 'squaredeuclidean');

        % Log-Gaussian kernel
        log_kn = -0.5 * D_sq;

        % Numerical stabilization
        max_log = max(log_kn, [], 2);
        log_sum = log(sum(exp(log_kn - max_log), 2)) + max_log;

        % Density estimation
        f_batch = exp(log_norm_const + log_sum - log(n));

        % Store results
        f(start_idx:end_idx) = gather(f_batch); % For GPU compatibility
    end
end