function fit_s(self)
    % increase accuracy of this fit
    % fmincon
    % spectrum
    self.order = '1';
    self.iiter = 0;
    lb = self.lb;
    ub = self.ub;
    x0 = self.x;

    nonlcon = [];
    [x] = nonlinear_fit(self, @(x) ObjFunBiXiAlpha(self, x), x0, lb, ub, nonlcon);
    [self.ks, self.kh, self.para_xi, self.para_alpha] = vecpara2para(x, self.order, self.para_xi, self.para_alpha, self.para_bs, self.model, self.peak_relation);

end