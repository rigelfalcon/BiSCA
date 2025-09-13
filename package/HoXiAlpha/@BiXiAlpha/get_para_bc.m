function [beta_bc, gamma_bc] = get_para_bc(self)
shat_xi = self.get_shat_component('xi');

f = self.f;
[~, idx_mu_xi] = max(shat_xi, [], 1);
mu_alpha = get_harmonic_mu_sigma(self.para_alpha.kernel.para.mu, 1, self.kh, self.df);
idx_mu_alpha = find_closest(f, mu_alpha, 1);

idx = get_ndgrid(repmat([idx_mu_xi(:); idx_mu_alpha(:)], 1, 2));
bchat = self.bchat;
ind = sub2ind(size(bchat), idx{1}, idx{2});
bchat = self.bchat(ind);
beta_bc = real(bchat);
gamma_bc = imag(bchat);
beta_bc(isinf(beta_bc)) = NaN;
gamma_bc(isinf(gamma_bc)) = NaN;

end
