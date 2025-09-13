function init_alpha(self)

    %% get init value
    init_para_alpha(self);
    bound_alpha(self);

    if strcmpi(self.method_alpha, 'findmax') && (strcmpi(self.peak_relation, 'free')) %isempty(self.bs) &&
        % findmax: locate peak positions (mu) freely from the smoothed spectrum.
        % Total peak count K = ks + kh is preserved; only the ks/kh split and
        % mu positions are determined by the data. K changes only if
        % find_largest_k_peaks cannot locate enough peaks in the spectrum.
        h = 0.05; % 0.03
        [s] = NwSmooth(self.f, mean(self.s, 2), h, self.f);
        min_distance = 2; % Minimum distance between peaks in frequency units

        [mu, h, baseline, s_detrended] = find_largest_k_peaks(self.f, s, self.k, min(self.f) + self.para_alpha.mu_df_lb, max(self.f) - self.para_alpha.mu_df_ub, min_distance);
        h(h < 0) = 0;

        % Find the index of the peak with the highest height
        [~, idx_max_h] = max(h);

        % Frequency of the peak with the highest height
        f_max_h = mu(idx_max_h);

        if ~self.para_alpha.no_zero_peak
            % Only add a DC "peak" if 0 is actually in the frequency grid.
            % Otherwise, mu will contain 0 but h computed via ismember(self.f,mu)
            % will not, causing dimension mismatches downstream.
            if ismember(0, self.f)
                mu = [0; mu];
                h = [0; h];
            end
            % Reclassify peaks as sub-/harmonic relative to the dominant peak.
            % This only redistributes the ks/kh split; total K = ks + kh is
            % preserved (equals the original requested K when all peaks are found).
            self.ks = sum(mu < f_max_h) - 1;
            self.ksmax = self.ks;
        else
            % Count the number of peaks with frequency smaller than f_max_h
            self.ks = sum(mu < f_max_h);
            self.ksmax = self.ks;
        end

        % Count the number of peaks with frequency larger than f_max_h
        self.kh = sum(mu >= f_max_h);
        self.khmax = self.kh;
        init_para_alpha(self);
        bound_alpha(self);
        self.para_alpha.kernel.para.mu = mu;
        shat_xi = pred_s_xi(self.f, self.para_xi);
        s_alpha = self.s - shat_xi;
        s_alpha_min = min(s_alpha);
        s_alpha = s_alpha - s_alpha_min;

        % Logging: mu must land on self.f grid, otherwise h will be shorter than mu.
        tf_mu_on_f = ismember(mu, self.f);
        if any(~tf_mu_on_f)
            mu_off = mu(~tf_mu_on_f);
            nshow = min(numel(mu_off), 8);
            warning('[BiXiAlpha:init_alpha] %d/%d mu not on f grid. Example(s): %s', ...
                numel(mu_off), numel(mu), mat2str(mu_off(1:nshow)));
        end
        tf_f_at_mu = ismember(self.f, mu);
        if sum(tf_f_at_mu) ~= numel(mu)
            warning('[BiXiAlpha:init_alpha] size mismatch risk: numel(mu)=%d but numel(h)=%d (mu on-grid=%d)', ...
                numel(mu), sum(tf_f_at_mu), sum(tf_mu_on_f));
        end

        self.para_alpha.kernel.para.h = s_alpha(tf_f_at_mu);

        if self.verbose
            figure(191);
            clf;
            subplot(2, 1, 1);
            hold on;
            plot(self.f, [self.s, shat_xi, s_alpha]);
            plot(mu, [self.para_alpha.kernel.para.h], 'ro', 'MarkerSize', 8, 'LineWidth', 2);
            legend(['s', 'shat_xi', 's_alpha']);
            figure(191);
            subplot(2, 1, 2);
            plot(self.f, [self.s, s, baseline, s_detrended]);
            hold on;
            baseline_at_mu = interp1(self.f, baseline, self.para_alpha.kernel.para.mu);
            plot(self.para_alpha.kernel.para.mu, h + baseline_at_mu, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
            legend(['s', 's_smooth', 'baseline', 's-detrended']);
        end

        self.para_alpha.mu_init = mu;
        self.para_alpha.sigma_init = 0.8;
        self.fx = [];
        self.fy = [];
        self.fxy = [];
        self.fxfy = [];
        self.para_bs.idx_fxfy = [];
    elseif strcmpi(self.method_alpha, 'findmax') && strcmpi(self.peak_relation, 'harmonic')
        h = 0.6;
        [s] = NwSmooth(self.f, mean(self.s, 2), h, self.f);
        range = [self.para_alpha.kernel.lb.mu(1), self.para_alpha.kernel.ub.mu(2)];

        [~, self.para_alpha.kernel.para.mu] = findmax_spectrum(log10(s - pred_s_xi(self.f, self.para_xi)), self.f, range); %
        self.para_alpha.mu_init = self.para_alpha.kernel.para.mu;
        self.para_alpha.sigma_init = 0.8;

        s_alpha = exp(mean(log(s), 2)) - pred_s_xi(self.f, self.para_xi);
        s_alpha_min = min(s_alpha);
        s_alpha = s_alpha - s_alpha_min;
        opt.maxopt = 'raw';
        self.para_alpha.kernel.para.h = find_max_around_xt(self.f, s_alpha, get_harmonic_mu_sigma(self.para_alpha.kernel.para.mu, 1, self.ks, self.kh, self.para_alpha.no_zero_peak)', min(self.width_candidate, self.para_alpha.kernel.para.mu / 2), opt) + s_alpha_min;
        self.para_alpha.kernel.para.h = self.para_alpha.kernel.para.h(:);
        self.para_alpha.kernel.para.h(1:self.k0 + self.ks) = 1e-1 * self.para_alpha.kernel.para.h(1:self.k0 + self.ks);
        self.para_alpha.isubharmonic = false;

        [self.para_alpha.kernel.para.mu] = get_harmonic_mu_sigma(self.para_alpha.kernel.para.mu, 1, self.ks, self.kh, self.para_alpha.no_zero_peak)';

        self.para_alpha.kernel.para.h(isnan(self.para_alpha.kernel.para.h) | self.para_alpha.kernel.para.h < 0) = 0;

        fit_alpha(self);

    elseif strcmpi(self.method_alpha, 'maxbs')

        if isempty(self.mufix)
            h = 0.6;
            [s] = NwSmooth(self.f, mean(self.s, 2), h, self.f);

            h = 0.6;
            bs_diag = abs(diag(mean(self.bs, 3)));
            bs_diag(isinf(abs(bs_diag))) = nan;
            bs_diag = NwSmooth(self.f, bs_diag, h, self.f);

            h = 0.6;
            bc_diag = abs(diag(mean(self.bs, 3)));
            bc_diag(isinf(abs(bc_diag))) = nan;
            bc_diag = NwSmooth(self.f, bc_diag, h, self.f);

            idx = find(self.f > range(1));
            bs_diag = bs_diag(idx);
            [maxbs, idx_max] = max(bs_diag);

            if idx_max == 1
                % might be the chaotic signal
                bc_diag = bc_diag(idx);
                [maxbc, idx_max] = max(bc_diag);
            end

            self.para_alpha.kernel.para.mu = self.f(idx(idx_max));

        else
            self.para_alpha.kernel.para.mu = self.mufix;
        end

        s_alpha = exp(mean(log(s), 2)) - pred_s_xi(self.f, self.para_xi);
        s_alpha_min = min(s_alpha);
        s_alpha = s_alpha - s_alpha_min;
        opt.maxopt = 'raw';
        self.para_alpha.kernel.para.h = find_max_around_xt(self.f, s_alpha, get_harmonic_mu_sigma(self.para_alpha.kernel.para.mu, 1, self.ks, self.kh, self.para_alpha.no_zero_peak)', min(self.width_candidate, self.para_alpha.kernel.para.mu / 2), opt) + s_alpha_min;
        self.para_alpha.kernel.para.h = self.para_alpha.kernel.para.h(:);
        self.para_alpha.kernel.para.h(1:self.k0 + self.ks) = 1e-1 * self.para_alpha.kernel.para.h(1:self.k0 + self.ks);
        self.para_alpha.isubharmonic = false;

    end

    self.para_alpha.kernel.para.h(isnan(self.para_alpha.kernel.para.h) | self.para_alpha.kernel.para.h < 0) = 0;

end