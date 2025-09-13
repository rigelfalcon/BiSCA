function [y_pred, low_dim] = nar_lpr_forecast(nar, init_state, steps, varargin)
    % NAR_LPR_FORECAST - Forecast using NAR-LPR model
    % Inputs:
    %   nar - Trained model from nar_lpr
    %   init_state - Initial state vector (1xp)
    %   steps - Number of steps to forecast
    % Outputs:
    %   y_pred - Forecasted values (steps x 1)
    %   low_dim - Low-dimensional projections for visualization

    % Parse optional parameters
    pin = inputParser;
    addParameter(pin, 'vis', false, @islogical);
    addParameter(pin, 'vis_dim', 3, @isnumeric);
    addParameter(pin, 'vis_train', false, @islogical);
    addParameter(pin, 'vis_forecast', true, @islogical);

    addParameter(pin, 'nburning', 0, @isnumeric);
    addParameter(pin, 'sigma', 0, @isnumeric);
    addParameter(pin, 'innovation_type', 'gaussian', @isCharString);
    addParameter(pin, 'innovation', [], @isnumeric);

    parse(pin, varargin{:});
    vis = pin.Results.vis;
    vis_dim = pin.Results.vis_dim;
    vis_train = pin.Results.vis_train;
    vis_forecast = pin.Results.vis_forecast;

    nburning = pin.Results.nburning;
    sigma = pin.Results.sigma;
    innovation_type = pin.Results.innovation_type;
    innovation = pin.Results.innovation;

    % Initialize variables
    y_pred = zeros(steps, size(init_state, 2));
    current_state = gpuArray(init_state);
    [N, p] = size(nar.X);

    % Get low-dimensional projection (PCA) from training data
    [coeff, ~, ~, ~, explained] = pca(nar.X);
    for i = 1:steps + nburning
        % Nadaraya-Watson prediction
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

        if ~isempty(sigma) && sigma > 0
            switch innovation_type
                case {'gaussian'}
                    noise = (nar.residual_chol * randn(nar.p, 1))';
                    y_next = y_next + sigma .* noise;
                case {'state_denpend'}
                    y_next = y_next + 0.6 * sigma .* innovation(i, :) + 0.4 * sigma .* innovation(i, :) .* mean(abs(y_next));
                case {'scaled_input'}
                    y_next = y_next + sigma .* innovation(i, :);
                case {'id_input'}
                    y_next = y_next + sigma .* (diag(diag(nar.residual_chol)) * innovation(i, :)')';
                case {'input'}
                    y_next = y_next + sigma .* (nar.residual_chol * innovation(i, :)')';
                otherwise
            end
        end

        % Store prediction
        y_pred(i, :) = y_next;

        % Update state (shift window)
        current_state = y_next;

    end
    y_pred = y_pred(nburning + 1:end, :);
    low_dim = y_pred * coeff(:, 1:vis_dim);

    % Visualization
    if vis
        hold on;

        % Project training data using correct dimensions
        train_proj = nar.Y * coeff(:, 1:vis_dim);

        % Choose plotting style based on visualization dimension
        if vis_dim == 2
            % 2D Visualization
            plot(train_proj(:, 1), train_proj(:, 2), 'b-', 'LineWidth', 1);
            h_pred = plot(low_dim(:, 1), low_dim(:, 2), 'r-', 'LineWidth', 2);
            scatter(low_dim(:, 1), low_dim(:, 2), 80, 'r', 'filled', ...
                'MarkerEdgeColor', 'k');

            xlabel('Principal Component 1');
            ylabel('Principal Component 2');

        else
            % 3D Visualization
            if vis_train
                if ~vis_forecast
                    plot3(train_proj(:, 1), train_proj(:, 2), train_proj(:, 3), ...
                        'b-'); %, 'LineWidth', 1
                else
                    plot3(train_proj(:, 1), train_proj(:, 2), train_proj(:, 3), ...
                        '.', 'Color', [0.5, 0.5, 0.5]); %, 'LineWidth', 1
                end

            end
            if vis_forecast
                h_pred = plot3(low_dim(:, 1), low_dim(:, 2), low_dim(:, 3), ...
                    'r-', 'LineWidth', 2);
                scatter3(low_dim(:, 1), low_dim(:, 2), low_dim(:, 3), ...
                    30, 'r', 'filled', 'MarkerEdgeColor', 'k');
            end
            xlabel('PC1');
            ylabel('PC2');
            zlabel('PC3');
            view(-30, 25)
        end

        grid off;
        axis tight;
        set(gca, 'FontSize', 12);
    end

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