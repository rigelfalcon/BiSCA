function fit_bs(self)

    %% remove the duplicate impact of wide range bispectrum/bicoherence from xi-xi xi+alpha
    if strcmpi(self.hos_components, 'xi+alpha')
        step = {'boundary'; 'diag'; 'offdiag'};
        self = fit_beta_gamma(self, step);
    elseif strcmpi(self.hos_components, 'alpha')
        step = {'diag'; 'offdiag'};
        self = fit_beta_gamma(self, step);
    elseif strcmpi(self.hos_components, 'xi')
        step = {'diag'};
        self = fit_beta_gamma(self, step);
    else
        error('no such method');
    end

    %% deal with nans
    self.para_bs.kernel.para.h(isnan(self.para_bs.kernel.para.h)) = 0;

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function self = fit_beta_gamma(self, step)
    % with xi and alpha
    % 1. get the beta and gamma of the xi-xi and xi+alpha from the boundary elements of the bispectrum first
    % 2. get the beta and gamma of the alpha-alpha from the diagonal elements of the bispectrum
    % 3. get the beta and gamma of the alpha-alpha from the off-diagonal elements of the bispectrum
    %
    self.order = '2';

    record.shat = self.shat;
    record.sxyhat = get_sxy(record.shat, self.f);
    record.sxsysxyhat = record.shat .* record.shat.' .* record.sxyhat;

    if any(ismember(step, 'boundary'))

        record.idx_mu = find_closest(self.f, 0);
        record.bs_boundary = self.bs(record.idx_mu, :)'; % boundary

        record.bsvar_boundary = self.bsvar(record.idx_mu, :)';
        record.bshatvar_boundary = self.bshatvar(record.idx_mu, :)';
        x0 = self.x;
        x0(isnan(x0)) = 0;

        self.region_bs = 'boundary';
        lb = self.lb;
        ub = self.ub;

        LostObjFunBoundary = @(x) ObjFunBiXiAlpha(self, x, record);
        [x] = nonlinear_fit(self, LostObjFunBoundary, x0, lb, ub);
        [~, ~, ~, ~, self.para_bs] = vecpara2para(x, self.order, self.para_xi, self.para_alpha, self.para_bs, self.model);

    end

    % fit the values from the diagonal
    if any(ismember(step, 'diag'))
        record.idx_diag = reshape((1:self.Nf) + (0:self.Nf - 1) * self.Nf, [], 1);
        record.bs_diag = self.bs(record.idx_diag); % diag

        record.bsvar_diag = self.bsvar(record.idx_diag);
        record.bshatvar_diag = self.bshatvar(record.idx_diag);
        x0 = self.x;
        x0(isnan(x0)) = 0;

        self.region_bs = 'diag';
        lb = self.lb;
        ub = self.ub;
        LostObjFunDiag = @(x) ObjFunBiXiAlpha(self, x, record);
        [x] = nonlinear_fit(self, LostObjFunDiag, x0, lb, ub);
        [~, ~, ~, ~, self.para_bs] = vecpara2para(x, self.order, self.para_xi, self.para_alpha, self.para_bs, self.model);

    end

    % fit the values from the off diagonal
    if any(ismember(step, 'offdiag'))
        record.bs = self.bs; % whole

        record.bsvar = self.bsvar;
        record.bshatvar = self.bshatvar;
        x0 = self.x;
        x0(isnan(x0)) = 0;

        self.region_bs = 'offdiag';
        lb = self.lb;
        ub = self.ub;
        LostObjFunOffDiag = @(x) ObjFunBiXiAlpha(self, x, record);
        [x] = nonlinear_fit(self, LostObjFunOffDiag, x0, lb, ub);
        [~, ~, ~, ~, self.para_bs] = vecpara2para(x, self.order, self.para_xi, self.para_alpha, self.para_bs, self.model);
    end

end