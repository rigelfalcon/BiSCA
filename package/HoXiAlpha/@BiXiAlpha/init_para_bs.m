function init_para_bs(self)

    switch self.para_bs.type
        case 'tstudent'
            if strcmpi(self.peak_relation, 'harmonic_free')
                if ~isfield(self.para_bs, 'kernel') || isempty(self.para_bs.kernel) || ~strcmpi(self.para_bs.kernel.type, self.kernel_type) || isempty(self.para_bs.kernel.para) || size(self.para_bs.kernel.para.h, 1) ~= self.N1MAX
                    self.para_bs.kernel = Kernel(self.kernel_type);

                    self.para_bs.kernel.para.h = ones(self.N1MAX, self.N1MAX, 2);
                    self.para_bs.kernel.para.mu = [self.para_xi.kernel.para.mu; self.para_alpha.kernel.para.mu(:)];
                    self.para_bs.kernel.para.sigma = [self.para_xi.kernel.para.sigma, 1; ...
                        self.para_alpha.kernel.para.sigma, ones(size(self.para_alpha.kernel.para.sigma))];
                    self.para_bs.kernel.para.nu = [self.para_xi.kernel.para.nu; self.para_alpha.kernel.para.nu];
                    self.para_bs.kernel.para.d = [self.para_xi.kernel.para.d; self.para_alpha.kernel.para.d];
                    self.para_bs.kernel.para.b = [self.para_xi.kernel.para.b; self.para_alpha.kernel.para.b];
                    self.para_bs.kernel.info.ks = self.ks;
                    self.para_bs.kernel.info.kh = self.kh;
                    self.para_bs.kernel.info.idx_sigma = 2 * ones(self.N1MAX, self.N1MAX);
                    self.para_bs.kernel.info.no_zero_peak = self.para_alpha.no_zero_peak;
                    self.para_bs.kernel.info.df = self.df;
                    self.para_bs.kernel.info.taperinfo = self.taperinfo;

                    if self.para_bs.no_xi_alpha
                        self.para_bs.kernel.info.idx_sigma(1, :) = 0;
                        self.para_bs.kernel.info.idx_sigma(:, 1) = 0;
                        self.para_bs.kernel.info.idx_sigma(1, 1) = 1;
                    else
                        self.para_bs.kernel.info.idx_sigma(1, :) = 3;
                        self.para_bs.kernel.info.idx_sigma(:, 1) = 4;
                        self.para_bs.kernel.info.idx_sigma(1, 1) = 1;
                    end

                end
            elseif strcmpi(self.peak_relation, 'free')
                self.para_bs.kernel = Kernel(self.kernel_type);

                self.para_bs.kernel.para.h = ones(self.N1MAX, self.N1MAX, 2);
                mu1 = [self.para_xi.kernel.para.mu; self.para_alpha.kernel.para.mu];
                self.para_bs.kernel.para.mu = get_ndgrid([mu1, mu1], 'list');

                kmax_alpha = numel(self.para_alpha.kernel.para.h);
                alpha_sigma = self.para_alpha.kernel.para.sigma;
                if isscalar(alpha_sigma), alpha_sigma = repmat(alpha_sigma, kmax_alpha, 1); end
                alpha_nu = self.para_alpha.kernel.para.nu;
                if isscalar(alpha_nu), alpha_nu = repmat(alpha_nu, kmax_alpha, 1); end
                alpha_d = self.para_alpha.kernel.para.d;
                if isscalar(alpha_d), alpha_d = repmat(alpha_d, kmax_alpha, 1); end
                alpha_b = self.para_alpha.kernel.para.b;
                if isscalar(alpha_b), alpha_b = repmat(alpha_b, kmax_alpha, 1); end

                sigma1 = [self.para_xi.kernel.para.sigma; alpha_sigma];
                self.para_bs.kernel.para.sigma = get_ndgrid([sigma1, sigma1], 'list');
                nu1 = [self.para_xi.kernel.para.nu; alpha_nu];
                self.para_bs.kernel.para.nu = get_ndgrid([nu1, nu1], 'list');
                d1 = [self.para_xi.kernel.para.d; alpha_d];
                self.para_bs.kernel.para.d = get_ndgrid([d1, d1], 'list');
                b1 = [self.para_xi.kernel.para.b; alpha_b];
                self.para_bs.kernel.para.b = get_ndgrid([b1, b1], 'list');

                self.para_bs.kernel.info.ks = self.ks;
                self.para_bs.kernel.info.kh = self.kh;
                self.para_bs.kernel.info.idx_sigma = 2 * ones(self.N1MAX, self.N1MAX);
                self.para_bs.kernel.info.no_zero_peak = self.para_alpha.no_zero_peak;
                self.para_bs.kernel.info.df = self.df;
                self.para_bs.kernel.info.taperinfo = self.taperinfo;
            elseif strcmpi(self.peak_relation, 'harmonic')
                self.para_bs.kernel = Kernel(self.kernel_type);

                self.para_bs.kernel.para.h = ones(self.N1MAX, self.N1MAX, 2);
                mu1 = [self.para_xi.kernel.para.mu; self.para_alpha.kernel.para.mu];
                self.para_bs.kernel.para.mu = get_ndgrid([mu1, mu1], 'list');

                kmax_alpha = numel(self.para_alpha.kernel.para.h);
                alpha_sigma = self.para_alpha.kernel.para.sigma;
                if isscalar(alpha_sigma), alpha_sigma = repmat(alpha_sigma, kmax_alpha, 1); end
                alpha_nu = self.para_alpha.kernel.para.nu;
                if isscalar(alpha_nu), alpha_nu = repmat(alpha_nu, kmax_alpha, 1); end
                alpha_d = self.para_alpha.kernel.para.d;
                if isscalar(alpha_d), alpha_d = repmat(alpha_d, kmax_alpha, 1); end
                alpha_b = self.para_alpha.kernel.para.b;
                if isscalar(alpha_b), alpha_b = repmat(alpha_b, kmax_alpha, 1); end

                sigma1 = [self.para_xi.kernel.para.sigma; alpha_sigma];
                self.para_bs.kernel.para.sigma = get_ndgrid([sigma1, sigma1], 'list');
                nu1 = [self.para_xi.kernel.para.nu; alpha_nu];
                self.para_bs.kernel.para.nu = get_ndgrid([nu1, nu1], 'list');
                d1 = [self.para_xi.kernel.para.d; alpha_d];
                self.para_bs.kernel.para.d = get_ndgrid([d1, d1], 'list');
                b1 = [self.para_xi.kernel.para.b; alpha_b];
                self.para_bs.kernel.para.b = get_ndgrid([b1, b1], 'list');

                self.para_bs.kernel.info.ks = self.ks;
                self.para_bs.kernel.info.kh = self.kh;
                self.para_bs.kernel.info.idx_sigma = 2 * ones(self.N1MAX, self.N1MAX);
                self.para_bs.kernel.info.no_zero_peak = self.para_alpha.no_zero_peak;
                self.para_bs.kernel.info.df = self.df;
                self.para_bs.kernel.info.taperinfo = self.taperinfo;
            end

        otherwise
            error('no such method')
    end

end