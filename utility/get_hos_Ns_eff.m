function [Ns_eff, ar_model] = get_hos_Ns_eff(x, hos_params, M, use_cache, use_kde)
    global NS_EFF_CACHE %#ok<*GVMIS>

    % Extract parameters
    window = hos_params.window;
    nfft = hos_params.nfft;
    overlap = hos_params.overlap;
    NW = hos_params.NW;
    N = length(x);
    if nargin < 4 || isempty(use_cache)
        use_cache = true;
    end
    if nargin < 5 || isempty(use_kde)
        use_kde = false;
    end
    % Check cache for matching parameters
    if use_cache && ~isempty(NS_EFF_CACHE)
        match_idx = find(...
            NS_EFF_CACHE.window == window & ...
            NS_EFF_CACHE.nfft == nfft & ...
            NS_EFF_CACHE.overlap == overlap & ...
            NS_EFF_CACHE.NW == NW & ...
            NS_EFF_CACHE.N == N, 1);
    else
        match_idx = [];
    end
    if ~isempty(match_idx)
        Ns_eff = NS_EFF_CACHE.Ns_eff(match_idx);
        ar_model = NS_EFF_CACHE.ar_model(match_idx);
        fprintf('Retrieved Ns_eff = %.2f for NW = %.2f from cache\n', Ns_eff, NW);
        return;
    end

    %%
    ar_model = fit_best_ar(gather(x), 1:100, 'levinson');
    nburning = 0;

    % *** MODIFICATION: Compute residuals for KDE ***
    if use_kde
        % Compute residuals: x(t) - predicted x(t) based on AR model
        innovations = randn(N + nburning, 1, 1, 'gpuArray');
        x_pred = gen_ar_ts(gpuArray(-permute(ar_model.A(2:end), [1, 3, 2])), ...
            zeros(N + nburning, 1, 1, 'gpuArray'), sqrt(ar_model.e) .* innovations);
        residuals = x - x_pred;
        % Ensure residuals are on CPU for KDE
        residuals = gather(residuals);
        % Fit KDE to residuals
        kde = fitdist(residuals, 'Kernel', 'Kernel', 'normal');
    end

    % Compute Ns_eff (no cache hit)
    if true %test_gpu([],[],true)
        if use_kde
            % *** MODIFICATION: Sample innovations from KDE on GPU ***
            innovations = gpuArray(random(kde, [N + nburning, 1, M]));
            x_sim = gen_ar_ts(gpuArray(-permute(ar_model.A(2:end), [1, 3, 2])), ...
                zeros(N + nburning, 1, M, 'gpuArray'), innovations);
        else
            innovations = randn(N + nburning, 1, M, 'gpuArray');
            x_sim = gen_ar_ts(gpuArray(-permute(ar_model.A(2:end), [1, 3, 2])), ...
                zeros(N + nburning, 1, M, 'gpuArray'), sqrt(ar_model.e) .* innovations);
        end

    else
        if use_kde
            % *** MODIFICATION: Sample innovations from KDE on CPU ***
            innovations = random(kde, [N + nburning, 1, M]);
        else
            innovations = randn(N + nburning, 1, M);
        end
        x_sim = gen_ar_ts(-permute(ar_model.A(2:end), [1, 3, 2]), ...
            zeros(N + nburning, 1, M), sqrt(ar_model.e) .* innovations);
    end
    x_sim = x_sim(nburning + 1:end, :, :, :); % Remove burning period, shape: (N, M)
    x_sim = permute(x_sim, [1, 2, 4, 3]);

    hos = HOS(x_sim, hos_params.Fs);
    hos.frange = hos_params.frange;
    hos.window = window;
    hos.nfft = nfft;
    hos.overlap = overlap;
    hos.NW = NW;
    if isfield(hos_params, 'num_random_segments')
        hos.num_random_segments = hos_params.num_random_segments;
    end
    hos.tf_return_seg = false;
    hos.tf_return_res = false;
    hos.tf_return_var = false;
    hos.normalization = 'haubrich';
    hos.method = 'pmtm';
    hos.bc;

    bc = hos.bc;
    valid_bc = bc(~isnan(bc));
    sum_bc2 = sum(abs(valid_bc).^2);
    count_bc = numel(valid_bc);
    mean_bc2 = sum_bc2 / count_bc;
    Ns_eff = 1 / mean_bc2;

    % Append to cache
    Ns = hos.Ns;
    new_row = table(window, nfft, overlap, NW, N, Ns_eff, Ns, ar_model);
    if use_cache
        NS_EFF_CACHE = [NS_EFF_CACHE; new_row];
    end
    fprintf('Computed and cached Ns_eff = %.2f for NW = %.2f, Ns = %.2f\n', Ns_eff, NW, hos.Ns);
end