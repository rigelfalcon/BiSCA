function nar_pruned = nar_lpr_screen_equilibrium_keys(nar, xq_equilibrium, varargin)
    % SCREEN_EQUILIBRIUM_KEYS - Prune NAR-LPR model based on equilibrium contributions
    % Inputs:
    %   nar - Trained NAR-LPR model from nar_lpr
    %   xq_equilibrium - Equilibrium query points (nq x p)
    % Optional Parameters:
    %   'threshold' - Minimum weight sum to retain a data point (default: 1e-10)
    % Output:
    %   nar_pruned - Pruned NAR model with irrelevant keys removed

    % Parse optional parameters
    pin = inputParser;
    addParameter(pin, 'threshold', 1e-10, @isnumeric);
    parse(pin, varargin{:});
    threshold = pin.Results.threshold;

    % Compute weights for all equilibrium points
    if ~nar.isadaptive
        % Fixed bandwidth case
        [~, L] = NwSmoothLogGaussian(nar.X, nar.Y, nar.hmin, gpuArray(xq_equilibrium));
    else
        % Adaptive case
        [~, ~, ~, L] = NwSmoothMassCenterBatch(nar.X, nar.Y, nar.kmin, xq_equilibrium, nar.h_i);

    end

    % Sum weights across all equilibrium points
    sum_weights = sum(gather(L), 1);

    % Identify relevant data points
    relevant_indices = sum_weights > threshold;

    % Create pruned model
    nar_pruned = nar;
    nar_pruned.X = nar.X(relevant_indices, :);
    nar_pruned.Y = nar.Y(relevant_indices, :);
    if nar.isadaptive
        nar_pruned.h_i = nar.h_i(relevant_indices, :);
        %% then use h_i and kmin will give different result.
    end

end