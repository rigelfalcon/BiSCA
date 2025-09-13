function [yq, selected_k, h_i, gcv_scores] = NwSmoothMassCenterGCV(x, y, kList, xq)
    % NwSmoothMassCenterGCV - Adaptive bandwidth kernel smoother with GCV for k selection
    %
    % Inputs:
    %   x      - Data points (n x dx matrix)
    %   y      - Response values (n x 1 vector)
    %   kList  - List of k values to evaluate for GCV (e.g., [1, 2, 5, 10, 20, 50, 100])
    %   xq     - Query points (nq x dx matrix, optional, default: x)
    %
    % Outputs:
    %   yq          - Smoothed values at query points using the selected k
    %   selected_k  - The k value that minimizes the GCV score
    %   gcv_scores  - GCV scores for each k in kList (optional output)

    % Handle default inputs
    if nargin < 4 || isempty(xq)
        xq = x;
    end

    [n, dx] = size(x); % Number of data points
    dk = length(kList); % Number of k values to evaluate

    % Precompute distances for all k values in kList
    max_k = max(kList); % Find the maximum k in kList
    [~, dist] = knnsearch(gather(x), gather(x), 'K', max_k + 1); % Compute distances on CPU for stability

    % Initialize GCV scores array
    gcv_scores = zeros(dk, 1);
    trace_avg = zeros(dk, 1);

    % Loop over each k in kList to compute GCV scores
    for i = 1:dk
        k = kList(i);

        % Extract precomputed distances for the current k
        h_i = dist(:, k + 1) ./ sqrt(dx); % Bandwidths for the current k
        h_i = max(h_i, min(h_i(h_i > 0)));

        % Robustness: If x is on the GPU, ensure h_i is also on the GPU before calling the batch function.
        if isgpuarray(x)
            h_i = gpuArray(h_i);
        end

        % Compute smoothed estimates and weight matrix at training points
        [yhat, ~, diagL] = NwSmoothMassCenterBatch(x, y, k, x, h_i);

        % Calculate Mean Squared Error (MSE)
        mse = mean((y - yhat).^2, "all");

        trace_avg(i) = mean(diagL);
        % Calculate GCV score: MSE / (1 - trace_avg)^2
        gcv_scores(i) = mse / (1 - trace_avg(i))^2;
    end

    % Find the k that minimizes the GCV score
    [~, idmin] = min(gcv_scores);
    selected_k = kList(idmin);
    h_i = sqrt((dist(:, selected_k + 1).^2) ./ dx); % Bandwidths for the selected k

    % Robustness: Ensure final h_i is on the correct device before the final call.
    if isgpuarray(x)
        h_i = gpuArray(h_i);
    end

    % Compute final smoothed estimates at query points using the selected k
    [yq] = NwSmoothMassCenterBatch(x, y, selected_k, xq, h_i);

    % Display results for user feedback
    disp(['Selected k: ', num2str(selected_k)]);
    disp(['Minimum GCV: ', num2str(gcv_scores(idmin))]);

end