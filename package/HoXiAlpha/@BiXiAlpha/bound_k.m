function [lb, ub] = bound_k(self)

    if strcmpi(self.ktype, 'fix')
        lb = [self.ks; self.kh];
        ub = [self.ks; self.kh];
    elseif strcmpi(self.ktype, 'free')
        lb = [0; min(ceil(max(self.f) / self.para_alpha.kernel.para.mu), self.khmax)];
        ub = [self.ksmax; self.khmax];
    elseif strcmpi(self.ktype, 'max')
        lb = [self.ksmax; self.khmax];
        ub = [self.ksmax; self.khmax];
    end
end