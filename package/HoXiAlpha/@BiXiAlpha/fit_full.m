function fit_full(self)

    % avoid redundant
    record.shat = self.shat;
    record.bs = self.bs;
    record.bc = self.bc;

    switch self.loss_type
        case {'bicoherence'}
            record.bsvar = self.bsdenom;
            record.bshatvar = self.bshatdenom;
        case {'HeteroEstVar', 'HeteroscedasticityEstVar'}
            record.bsvar = self.bsvar;
            record.bshatvar = self.bshatvar;
    end

    self.region_bs = 'all';
    lb = self.lb;
    ub = self.ub;
    x0 = self.x;

    [x] = nonlinear_fit(self, @(x) ObjFunBiXiAlpha(self, x, record), x0, lb, ub);
    [~, ~, ~, ~, self.para_bs] = vecpara2para(x, self.order, self.para_xi, self.para_alpha, self.para_bs, self.model, self.peak_relation);

    % spectrum and bispectrum
    self.order = '1+2';
    self.region_bs = 'all';
    lb = self.lb;
    ub = self.ub;
    x0 = self.x;

    [x] = nonlinear_fit(self, @(x) ObjFunBiXiAlpha(self, x, record), x0, lb, ub);
    [self.ks, self.kh, self.para_xi, self.para_alpha, self.para_bs] = vecpara2para(x, self.order, self.para_xi, self.para_alpha, self.para_bs, self.model, self.peak_relation);

    self.model.jacobian = [];

end