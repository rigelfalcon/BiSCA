function H = mvts2hankel(x, p)
    % MVTS2HANKEL Construct multichannel Hankel matrix with channel-first, lag-second ordering
    % Inputs:
    %   x: [Nt × Nc] matrix (Nt time points, Nc channels)
    %   p: embedding dimension (number of lags)
    % Output:
    %   H: [Nt - p + 1 × Nc * p] matrix with columns ordered as:
    %      [Channel1_Lag1, Channel2_Lag1, ..., ChannelNc_Lag1,
    %       Channel1_Lag2, Channel2_Lag2, ..., ChannelNc_Lagp]

    % Ensure input is column-oriented
    if isrow(x)
        x = x';
    end
    [Nt, Nc] = size(x);

    % Validate embedding dimension
    if p > Nt
        error('Embedding dimension p exceeds number of time points');
    end

    % Preallocate output matrix
    num_windows = Nt - p + 1;
    H = zeros(num_windows, Nc * p, 'like', x);

    % Build Hankel matrix with channel-first ordering
    for c = 1:Nc
        % Get Hankel matrix for current channel [p × num_windows]
        chan_hankel = hankel(x(1:p, c), x(p:end, c));

        % Transpose to [num_windows × p], then assign columns by lag
        chan_hankel = chan_hankel';

        for lag = 1:p
            % Calculate column index: (lag-1)*Nc + current_channel
            col_idx = (lag - 1) * Nc + c;
            H(:, col_idx) = chan_hankel(:, lag);
        end
    end
end