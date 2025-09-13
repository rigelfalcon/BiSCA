function armodel = fit_best_ar(y, Nlag, method)
    % FIT_BEST_AR Fits AR models of different orders to data y and selects the best model based on AICc
    %
    %   armodel = FIT_BEST_AR(y, Nlag, method)
    %
    %   Inputs:
    %       y      - Time series data (vector)
    %       Nlag   - Vector of AR model orders to try
    %       method - (Optional) Method to estimate AR parameters. Options are:
    %                'arburg'   - Uses Burg's method
    %                'levinson' - Uses Levinson-Durbin recursion (default)
    %
    %   Output:
    %       armodel - Struct containing the best AR model (fields: A, e)
    %

    % Set default method to 'levinson' if not provided
    if nargin < 3
        method = 'levinson';
    end

    % Ensure y is a column vector
    if size(y, 2) > 1
        y = y(:);
    end

    N = length(y);

    % Ensure Nlag is a column vector and sort in ascending order
    Nlag = Nlag(:);
    Nlag = sort(Nlag);

    % Convert y to double if it's single precision
    if isa(y, 'single')
        y = double(y);
    end

    % Initialize models cell array
    models = cell(length(Nlag), 1);

    % Start timing
    tic

    % Select method to compute AR parameters
    switch lower(method)
        case 'arburg'
            % Use Burg's method
            for i = 1:length(Nlag)
                p = Nlag(i);
                [A, E] = arburg(y, p);
                models{i}.A = A; % Normalized AR coefficients (leading 1)
                models{i}.e = E; % Estimated white noise variance
            end

        case 'levinson'
            % Use Levinson-Durbin recursion
            % Compute autocorrelation sequence once up to maximum lag
            maxLag = max(Nlag);
            r_full = xcorr(y, maxLag, 'biased');
            % Extract positive lags (from lag 0 to maxLag)
            r = r_full(maxLag + 1:end);

            for i = 1:length(Nlag)
                p = Nlag(i);
                % Apply Levinson-Durbin recursion
                [A, E] = levinson(r(1:p + 1), p);
                models{i}.A = A; % AR coefficients
                models{i}.e = E; % Prediction error (variance of white noise input)
            end

        otherwise
            error('Unknown method "%s". Options are "arburg" or "levinson".', method);
    end

    % End timing
    t = toc;
    fprintf('Fitting AR models took %f seconds\n', t);

    % Compute AICc for each model
    V = nan(length(Nlag), 1);

    for i = 1:length(Nlag)
        V(i) = araicc(models{i}, N);
    end

    % Select the model with minimum AICc
    [~, I] = min(V);
    armodel = models{I};
    fprintf('Choose order  %f \n', Nlag(I));

end

% AICc computation function
function v = araicc(model, N)
    % ARAICC Computes the corrected Akaike Information Criterion (AICc) for an AR model
    %
    %   v = ARAICC(model, N)
    %
    %   Inputs:
    %       model - Struct containing AR model parameters (fields: A, e)
    %       N     - Length of the data used to fit the model
    %
    %   Output:
    %       v - AICc value
    %

    p = length(model.A) - 1; % Model order
    k = p + 1; % Number of parameters (AR coefficients + variance)

    % Log-likelihood
    LL = -0.5 * N * log(model.e);

    % Akaike Information Criterion
    AIC = 2 * k - 2 * LL;

    % Corrected AIC
    v = AIC + (2 * k * (k + 1)) / (N - k - 1);

end