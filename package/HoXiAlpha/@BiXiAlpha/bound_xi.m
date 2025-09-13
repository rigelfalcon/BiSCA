function [lb, ub] = bound_xi(self)

    idx_xi_range = self.f < self.para_xi.f_max;

    if strcmpi(self.para_xi.type, 'tstudent')
        self.para_xi.kernel.lb.h = max(mean(self.s(idx_xi_range, :), 2)) / 100;
        self.para_xi.kernel.lb.mu = 0;
        self.para_xi.kernel.lb.sigma = 0.5; %1 1e-5 15 1;
        self.para_xi.kernel.lb.nu = 1e-4; %0.1

        if ~self.para_xi.use_xid && isfield(self.para_alpha, 'kernel') && isfield(self.para_alpha.kernel, 'lb') && isfield(self.para_alpha.kernel.lb, 'd')
            self.para_xi.kernel.lb.d = self.para_alpha.kernel.lb.d;
        else
            self.para_xi.kernel.lb.d = 1;
        end

        self.para_xi.kernel.lb.b = 1; %0
        self.para_xi.kernel.ub.h = 1.5 * max(self.s, [], 'all');
        self.para_xi.kernel.ub.mu = 0;
        self.para_xi.kernel.ub.sigma = 1e4;
        self.para_xi.kernel.ub.nu = 1e2; %5

        if ~self.para_xi.use_xid && isfield(self.para_alpha, 'kernel') && isfield(self.para_alpha.kernel, 'ub') && isfield(self.para_alpha.kernel.ub, 'd')
            self.para_xi.kernel.ub.d = self.para_alpha.kernel.ub.d;
        else
            % d controls flat-top width: larger d -> narrower flat-top -> more AR(1)-like decay
            % Relaxed from 10 to 1e3 to allow better fit to low-pass spectra
            self.para_xi.kernel.ub.d = 1e3;
        end

        self.para_xi.kernel.ub.b = 1; %1e5

        if ~isempty(self.para_xi.kernel.fix)
            self.para_xi.kernel.vlb(~isnan(self.para_xi.kernel.vfix)) = self.para_xi.kernel.vfix(~isnan(self.para_xi.kernel.vfix));
            self.para_xi.kernel.vub(~isnan(self.para_xi.kernel.vfix)) = self.para_xi.kernel.vfix(~isnan(self.para_xi.kernel.vfix));
        end

        lb = self.para_xi.kernel.vlb(:);
        ub = self.para_xi.kernel.vub(:);

    elseif strcmpi(self.para_xi.type, 'tstudent+peak')
        lb = [max(self.s(idx_xi_range, :)) / 10, 4, 0, 0, ... %xi: A B nu U
            0, min(self.df, 0.5), 0]; %xi: mu sigma alpha
        ub = [2 * max(self.s, [], 'all'), 100, 10, 0, ... %xi: A B nu U
            0, 4, max(self.s, [], 'all')]; %xi: mu sigma alpha
    elseif strcmpi(self.para_xi.type, 'gaussian')
        self.para_xi.kernel.lb.alpha = max(self.s(idx_xi_range, :)) / 10;
        self.para_xi.kernel.lb.mu = 0;
        self.para_xi.kernel.lb.sigma = 4;
        self.para_xi.kernel.ub.alpha = 2 * max(self.s, [], 'all');
        self.para_xi.kernel.ub.mu = 0;
        self.para_xi.kernel.ub.sigma = 100;

        lb = self.para_xi.kernel.vlb(:)';
        ub = self.para_xi.kernel.vub(:)';
    elseif strcmpi(self.para_xi.type, 'lorentzian')
        self.para_xi.kernel.lb.h = max(mean(self.s(idx_xi_range, :), 2)) / 100;
        self.para_xi.kernel.lb.mu = 0;
        self.para_xi.kernel.lb.k = 0; %1 1e-5 15 1;
        self.para_xi.kernel.lb.chi = 1e-1; %0.1

        self.para_xi.kernel.ub.h = 1.5 * max(self.s, [], 'all');
        self.para_xi.kernel.ub.mu = 0;
        self.para_xi.kernel.ub.k = 1e3;
        self.para_xi.kernel.ub.chi = 1e2; %5

        if ~isempty(self.para_xi.kernel.fix)
            self.para_xi.kernel.vlb(~isnan(self.para_xi.kernel.vfix)) = self.para_xi.kernel.vfix(~isnan(self.para_xi.kernel.vfix));
            self.para_xi.kernel.vub(~isnan(self.para_xi.kernel.vfix)) = self.para_xi.kernel.vfix(~isnan(self.para_xi.kernel.vfix));
        end

        lb = self.para_xi.kernel.vlb(:);
        ub = self.para_xi.kernel.vub(:);

    elseif strcmpi(self.para_xi.type, 'linear')
        lb = [-1, 0];
        ub = [0, max(self.s, [], 'all')];
    else
        error('no such self.para_xi.type')
    end

end