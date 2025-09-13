function x_resampled = resampleGridded(x_original, Fs_original, Fs_upsampled)
    % Determine the original time vector
    N_original = size(x_original, 1);
    t_original = ((0:N_original - 1).') / Fs_original;
    t_original_end = t_original(end);

    % Calculate the upsampled time vector
    step = 1 / Fs_upsampled;
    t_upsampled = (0:step:t_original_end).';

    % Ensure the upsampled time does not exceed t_original_end due to floating-point precision issues
    t_upsampled = t_upsampled(t_upsampled <= t_original_end);

    % Create the gridded interpolant with linear interpolation and no extrapolation (returns NaN)
    F = griddedInterpolant(t_original, x_original, 'cubic', 'none');

    % Perform interpolation on the upsampled time grid
    x_resampled = F(t_upsampled);
end