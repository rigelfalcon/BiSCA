function fit_xi(self)

    x0 = self.para_xi.kernel.vpara;

    [lb, ub] = bound_xi(self);

    x = nonlinear_fit(self, @(x) LostObjFun(x, self), x0, lb, ub);
    self.para_xi.kernel.vpara = x;
    self.trend = pred_s_xi(self.f, self.para_xi);
end
%% lost function for the xi process
function loss = LostObjFun(x, self)
    x = recover_fixed_para(x, [], [], self.para_fit.parainfo);
    para_xi = self.para_xi;

    para_xi.kernel.vpara = x;

    shat_xi = pred_s_xi(self.f, para_xi);
    self.iiter = self.iiter + 1;
    if self.verbose && mod(self.iiter, 100) == 1
        figure(31)
        plot(self.f, log10(abs([self.s, shat_xi])))
        drawnow
    end

    % Use f_max to focus xi fitting on low frequencies where xi dominates
    % This prevents xi from absorbing power from alpha peaks at higher frequencies
    f_max = para_xi.f_max;  % Default is 3 Hz from init_para_xi
    idx_low = self.f <= f_max;
    idx_high = self.f > f_max;
    
    % Primary loss: fit xi spectrum only in low-frequency range (f <= f_max)
    % where xi is expected to dominate the composite signal
    loss_low = calc_loss_s(self, self.loss_type, self.s(idx_low), shat_xi(idx_low), self.svar(idx_low), []);
    
    % Endpoint constraints
    diff_endpoint = [self.s(1) - shat_xi(1); self.s(end) - shat_xi(end)];
    
    % Trend residual penalty: xi should not exceed observed spectrum anywhere
    % This is applied to ALL frequencies to ensure xi stays below composite
    trend_residual = shat_xi - self.s;
    trend_residual(trend_residual < 0) = 0;
    
    % Additional penalty for high frequencies: discourage xi power above f_max
    % The xi process should decay above f_max, not form a "shoulder"
    % Use a moderate penalty (100) to allow some high-freq power but discourage excess
    high_freq_penalty = shat_xi(idx_high) ./ max(self.s(idx_high));  % Normalized
    
    penalty = [1e3 .* diff_endpoint(1); ... % start point
        2 * 1e3 .* diff_endpoint(2); ...    % end point
        1e4 .* trend_residual; ...          % xi <= s everywhere
        1e3 .* high_freq_penalty];          % discourage xi power at high freq
    loss = [loss_low; penalty];
    self.model.loss_history(self.iiter, :) = sum(loss(~isnan(loss)));
end
