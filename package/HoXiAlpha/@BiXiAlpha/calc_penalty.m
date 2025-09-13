function penalty = calc_penalty(self, para_xi, para_alpha, para_bs)

    if ~isempty(self.lambda) && ((isscalar(self.lambda) && self.lambda ~= 0) || length(self.lambda) > 1)

        switch self.regularization
            case {''}
                penalty = 0;
            case {'L2', 'l2'}

                penalty = self.lambda(1) * norm([para_xi.kernel.para.h(:); para_alpha.kernel.para.h(:)]) + self.lambda(2) * norm([para_bs.kernel.para.h(:)]);
            case {'L1', 'l1'}

                penalty = self.lambda(1) * norm([para_xi.kernel.para.h(:); para_alpha.kernel.para.h(:)], 1) + self.lambda(2) * norm([para_bs.kernel.para.h(:)], 1);

            case {'L2,L1', 'l2,l1'}

                penalty = self.lambda(1) * norm([para_xi.kernel.para.h(:); para_alpha.kernel.para.h(:)], 2) + self.lambda(2) * norm([para_bs.kernel.para.h(:)], 1);

            case {'L1-subharmonic'}

                if strcmpi(self.order, '2') || strcmpi(self.order, '1+2')

                    penalty = self.lambda .* norm([reshape(para_bs.kernel.para.h(2:self.ks + self.k0 + 1, :, :), [], 1); reshape(para_bs.kernel.para.h(:, 2:self.ks + self.k0 + 1, :), [], 1)], 1);

                else
                    penalty = 0;
                end
            case {'elastic_constraint_xi'}

                shat_xi = pred_s_xi(self.f, para_xi);
                diff_endpoint = [(self.s(1) - shat_xi(1)) ./ self.s(1); (self.s(end) - shat_xi(end)) ./ self.s(end)];
                trend_residual = shat_xi - self.s;
                trend_residual(trend_residual < 0) = 0;
                penalty = [self.lambda(1) .* diff_endpoint(1); ... %start
                    self.lambda(2) .* diff_endpoint(2); ... %end
                    self.lambda(3) .* trend_residual]; % allways below empirical s
            case {'elastic_constraint_mu_sigma_chi_xi'}
                diff_mu = self.para_alpha.kernel.para.mu(:) - self.para_alpha.mu_init(:);
                diff_sigma = para_alpha.kernel.para.sigma - self.para_alpha.sigma_init;
                diff_sigma(diff_sigma < 0) = 0;
                shat_xi = pred_s_xi(self.f, para_xi);
                diff_endpoint = [(self.s(1) - shat_xi(1)) ./ self.s(1); (self.s(end) - shat_xi(end)) ./ self.s(end)];
                trend_residual = shat_xi - self.s;
                trend_residual(trend_residual < 0) = 0;
                penalty = [self.lambda(1) .* diff_mu; ... %mu
                    self.lambda(2) .* diff_sigma; ... %sigma
                    self.lambda(3) .* diff_endpoint(1); ... %start
                    self.lambda(4) .* diff_endpoint(2); ... %end
                    self.lambda(5) .* trend_residual]; % allways below empirical s
            otherwise
                error('no such method')
        end

    else
        penalty = 0;
    end

end