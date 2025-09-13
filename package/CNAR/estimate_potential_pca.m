function result = estimate_potential_pca(nar, num_steps, grid_points, sigma, M)
    % ESTIMATE_POTENTIAL_PCA - Computes potential landscape and related data
    % Inputs:
    %   nar         - NAR model
    %   num_steps   - Number of simulation steps
    %   grid_points - Number of grid points per dimension
    %   sigma       - Noise magnitude
    %   M           - Number of parallel trials
    % Output:
    %   result      - Struct containing all computed data

    % 1. Simulate in original space
    x_high = simulate_high_dim(nar, num_steps, sigma, M);

    % Upsample trajectory
    Fs_original = 200;
    Fs_upsampled = 2000;
    x_high = reshape(resampleGridded(reshape(x_high, size(x_high, 1), []), Fs_original, Fs_upsampled), [], size(x_high, 2), size(x_high, 3));

    % 2. Perform PCA
    [coeff, ~, ~, ~, explained] = pca(nar.X);

    % 3. Project to PCA space
    x_low = project_to_low_dim(x_high, coeff);

    % 4. Create grid based on data range
    min_x = min(x_low(:, 1)) - 0.1;
    max_x = max(x_low(:, 1)) + 0.1;
    min_y = min(x_low(:, 2)) - 0.1;
    max_y = max(x_low(:, 2)) + 0.1;
    grid_vals = {linspace(min_x, max_x, grid_points), linspace(min_y, max_y, grid_points)};

    % 5. Compute optimal bandwidth using GCV
    X_low = nar.X * coeff(:, 1:2);
    Y_low = nar.Y * coeff(:, 1:2);
    hlist = get_hList(100, [(max(X_low, [], 'all') - min(X_low, [], 'all')) / 200, ...
        (max(X_low, [], 'all') - min(X_low, [], 'all')) / 1], @logspace);
    [~, gcv_low] = NwSmoothGCVOriginalScaleLogGaussian(X_low, Y_low, hlist, X_low);
    h = gcv_low.hmin / 2;

    % 6. Compute potential landscape
    [U, X1, X2] = compute_potential_from_low_dim(x_low, grid_vals, h);

    % 7. Package results into struct
    result = struct();
    result.X1 = X1;
    result.X2 = X2;
    result.U = U;
    result.x_low = reshape(x_low, size(x_high, 1), M, 2);
    result.coeff = coeff;
    result.explained = explained;
    result.grid_vals = grid_vals;
    result.h = h;
    result.num_steps = num_steps;
    result.sigma = sigma;
    result.M = M;
end

function x_high = simulate_high_dim(nar, num_steps, sigma, M, initial_states)
    % SIMULATE_HIGH_DIM - Vectorized high-dim stochastic simulations
    % Inputs:
    %   nar            - Original NAR model
    %   num_steps      - Steps per simulation
    %   sigma          - Noise magnitude
    %   M              - Number of parallel trials
    %   initial_states - Optional initial states (M x p matrix)
    % Output:
    %   x_high         - Simulated trajectories (num_steps x M x 1)

    nburning = 0;

    if nargin < 5
        % Random initial states from training data
        if M > size(nar.X, 1)
            idx = randi(size(nar.X, 1), M, 1);
        else
            idx = randperm(size(nar.X, 1), M);
        end
        initial_states = nar.X(idx, :);
    end

    p = size(nar.X, 2); % NAR model order
    x_high = zeros(num_steps + nburning, M, p);
    current_state = initial_states; % M x p

    % Vectorized simulation loop
    for t = 1:(num_steps + nburning)
        % Predict next value for all trials
        if ~nar.isadaptive
            [y_next] = NwSmoothLogGaussian(nar.X, nar.Y, nar.hmin, current_state);
        else
            [y_next] = NwSmoothMassCenterBatch(nar.X, nar.Y, nar.kmin, current_state, nar.h_i);
        end

        if isfield(nar, 'issym') && nar.issym
            if ~nar.isadaptive
                [y_next2] = NwSmoothLogGaussian(nar.X, nar.Y, nar.hmin, -current_state);
            else
                [y_next2] = NwSmoothMassCenterBatch(nar.X, nar.Y, nar.kmin, -current_state, nar.h_i);
            end
            % Enforce antisymmetry: f(-x) = -f(x)
            y_next = (y_next - y_next2) / 2;
        end

        % Add noise and store
        noise = (nar.residual_chol * randn(nar.p, 1))';
        y_next = y_next + sigma .* noise;

        x_high(t, :, :) = y_next;
        current_state = y_next;
    end
    x_high = x_high(nburning + 1:end, :, :);
end

function x_low = project_to_low_dim(x_high, coeff)
    % PROJECT_TO_LOW_DIM - Projects simulations to PCA space
    % Inputs:
    %   x_high - High-dim data (num_steps x M x 1)
    %   coeff  - PCA coefficients from reduced model
    % Output:
    %   x_low  - Low-dim projections (N x 2)

    % Reshape and project
    x_flat = reshape(x_high, [], size(x_high, 3)); % N x 1
    x_low = x_flat * coeff(:, 1:2); % N x 2
end

function [U, X1, X2] = compute_potential_from_low_dim(x_low, grid_vals, h)
    % COMPUTE_POTENTIAL_FROM_LOW_DIM - Calculates landscape from projections
    % Inputs:
    %   x_low     - Low-dim data (N x 2)
    %   grid_vals - Grid definition {x_grid, y_grid}
    % Outputs:
    %   U, X1, X2 - Potential landscape and meshgrid

    % Kernel density estimation
    [X1, X2] = meshgrid(grid_vals{:});
    [pdf] = kde_batch(x_low, [X1(:), X2(:)], h);

    % Potential calculation
    P_ss = reshape(pdf, size(X1));
    U = -log(P_ss + 1e-10);
    U = U - min(U(:)); % Zero-base potential
end

function [U, X1, X2] = compute_potential_from_low_dim_adptive(x_low, grid_vals)
    % COMPUTE_POTENTIAL_FROM_LOW_DIM - Calculates landscape from projections
    % Inputs:
    %   x_low     - Low-dim data (N x 2)
    %   grid_vals - Grid definition {x_grid, y_grid}
    % Outputs:
    %   U, X1, X2 - Potential landscape and meshgrid

    % Kernel density estimation
    [X1, X2] = meshgrid(grid_vals{:});
    klist = unique(ceil(linspace(1, 200, 100))');
    [pdf] = kde_masscenter_gcv(x_low, klist, [X1(:), X2(:)]);

    % Potential calculation
    P_ss = reshape(pdf, size(X1));
    U = -log(P_ss + 1e-10);
    U = U - min(U(:)); % Zero-base potential
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

    if nargout > 1
        dbg.s = sum(weights, 1);
        L = weights;
    end
end