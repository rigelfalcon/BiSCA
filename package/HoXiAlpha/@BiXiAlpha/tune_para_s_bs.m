function tune_para_s_bs(self)

    % only the diagonal and cross line of fundamental peak
    self.order = '1+2';

    record.idx_diag = reshape((1:self.Nf) + (0:self.Nf - 1) * self.Nf, [], 1);
    record.bs_diag = [self.bs(record.idx_diag)]; %self.bs(idx_mu, :)';
    record.bsvar_diag = [self.bsvar(record.idx_diag)]; %self.bshatvar(idx_mu, :)';
    record.bshatvar_diag = [self.bshatvar(record.idx_diag)]; %self.bshatvar(idx_mu, :)';

    % x0
    x0 = self.x;
    x0(isnan(x0)) = 0;

    self.region_bs = 'mu+diag';
    lb = self.lb;
    ub = self.ub;

    % fit
    LostObjFunTune = @(x) ObjFunBiXiAlpha(self, x, record);
    nonlcon = [];
    [x] = nonlinear_fit(self, LostObjFunTune, x0, lb, ub, nonlcon);
    [self.ks, self.kh, self.para_xi, self.para_alpha, self.para_bs] = vecpara2para(x, self.order, self.para_xi, self.para_alpha, self.para_bs, self.model);

    %% deal with nans
    self.para_bs.kernel.para.h(isnan(self.para_bs.kernel.para.h)) = 0;

end