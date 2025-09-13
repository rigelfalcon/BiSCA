function R2w = calc_r2w(x, x_hat, x_variance)
    %CALC_R2W Calculate Weighted R-squared for Observed and Predicted Data
    %
    %   R2w = calc_r2w(x, x_hat, x_variance) computes the weighted R-squared
    %   between the observed data `x` and predicted data `x_hat`, given the
    %   variances `x_variance`. The weights are defined as w = 1 ./ x_variance.
    %
    %   Inputs:
    %       x           - Observed data (numeric array, real)
    %       x_hat       - Predicted data (same size as x)
    %       x_variance  - Variance of each observation (same size as x)
    %
    %   Output:
    %       R2w         - Weighted R-squared (scalar)
    %
    %   Example Usage:
    %
    %       % --- Example 1: Real Data ---
    %       x = [10; 15; 20; 25; 30];
    %       x_hat = [12; 14; 19; 24; 28];
    %       x_variance = [2; 3; 2.5; 4; 3.5];
    %       R2w = calc_r2w(x, x_hat, x_variance);
    %       fprintf('Weighted R²: %.4f\n', R2w);
    %
    %       % --- Example 2: Matrix Data ---
    %       x = [1, 2, 3; 4, 5, 6];
    %       x_hat = [1.1, 1.9, 3.2; 3.8, 5.1, 6.3];
    %       x_variance = [0.1, 0.2, 0.3; 0.2, 0.1, 0.4];
    %       R2w = calc_r2w(x, x_hat, x_variance);
    %       fprintf('Weighted R²: %.4f\n', R2w);
    %
    %       % --- Example 3: Applying to Complex Data ---
    %       M = 50; N = 50;
    %       real_x = randn(M, N);
    %       imag_x = randn(M, N);
    %       x_complex = complex(real_x, imag_x); % Observed Bispectrum
    %
    %       real_x_hat = randn(M, N);
    %       imag_x_hat = randn(M, N);
    %       x_hat_complex = complex(real_x_hat, imag_x_hat); % Predicted Bispectrum
    %
    %       x_variance = rand(M, N) + 1; % Ensure positive variances
    %
    %       % Calculate Weighted R² for Real Part
    %       R2w_real = calc_r2w(real(x_complex), real(x_hat_complex), x_variance);
    %       fprintf('Weighted R² (Real Part): %.4f\n', R2w_real);
    %
    %       % Calculate Weighted R² for Imaginary Part
    %       R2w_imag = calc_r2w(imag(x_complex), imag(x_hat_complex), x_variance);
    %       fprintf('Weighted R² (Imaginary Part): %.4f\n', R2w_imag);
    %
    %       % Combined R² can be computed externally as needed
    %       R2w_combined = (R2w_real + R2w_imag) / 2;
    %       fprintf('Weighted R² (Combined): %.4f\n', R2w_combined);
    %

    %% 1. Input Validation
    % Ensure inputs are numeric arrays of the same size
    if ~isequal(size(x), size(x_hat), size(x_variance))
        error('Inputs x, x_hat, and x_variance must have the same dimensions.');
    end

    % Ensure all inputs are real-valued
    if ~isreal(x) || ~isreal(x_hat) || ~isreal(x_variance)
        error('All inputs x, x_hat, and x_variance must be real-valued.');
    end

    % Ensure x_variance contains positive values to avoid division by zero
    if any(x_variance(:) <= 0)
        error('All elements of x_variance must be positive.');
    end

    %% 2. Compute Weights
    w = 1 ./ x_variance; % Element-wise weights

    %% 3. Compute Weighted Mean of Observed Data
    weighted_mean_x = sum(w(:) .* x(:), "all", "omitnan") / sum(w(:), "all", "omitnan");

    %% 4. Compute Weighted Sum of Squared Errors (SSEw)
    residuals = x_hat - x;
    SSEw = sum(w(:) .* residuals(:).^2, "all", "omitnan");

    %% 5. Compute Weighted Total Sum of Squares (SSTw)
    SSTw = sum(w(:) .* (x(:) - weighted_mean_x).^2, "all", "omitnan");

    %% 6. Calculate Weighted R-squared
    R2w = 1 - (SSEw / SSTw);
end