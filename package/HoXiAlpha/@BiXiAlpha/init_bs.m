function init_bs(self)

    init_para_bs(self);
    bound_bs(self);
    bs = mean(self.bs, 3);
    bsre = real(bs);
    bsim = imag(bs);

    if strcmpi(self.peak_relation, 'harmonic_old')
        mu1 = [self.para_xi.kernel.para.mu; get_harmonic_mu_sigma(self.para_alpha.kernel.para.mu, 0, self.ks, self.kh, self.para_alpha.no_zero_peak)'];
        mu2 = get_ndgrid([mu1(:), mu1(:)], 'list');
    else
        mu1 = [self.para_xi.kernel.para.mu; self.para_alpha.kernel.para.mu(:)];
        mu2 = get_ndgrid([mu1(:), mu1(:)], 'list');
    end
    self.para_bs.kernel.para.h(:, :, 1) = reshape(find_max_around_xt(cat(3, self.fx, self.fy), bsre, mu2, self.width_candidate), self.N1, self.N1);
    self.para_bs.kernel.para.h(:, :, 2) = reshape(find_max_around_xt(cat(3, self.fx, self.fy), bsim, mu2, self.width_candidate), self.N1, self.N1);

    self.para_bs.kernel.para.h(isnan(self.para_bs.kernel.para.h)) = 0;

    idx_triu = logical(mat2triu(true(size(self.para_bs.kernel.para.h)), 1, false));
    self.para_bs.kernel.para.h(idx_triu) = 0;
end