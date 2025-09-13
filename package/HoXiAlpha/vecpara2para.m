% extract the parameter back from vectorized parameter
function [ks, kh, para_xi, para_alpha, para_bs] = vecpara2para(x, order, para_xi, para_alpha, para_bs, model, peak_relation)
    if strcmpi(peak_relation, 'harmonic_old')
        x = x(:);

        if nargin > 5 && ~isempty(model)
            return_stat = true;
        else
            return_stat = false;
        end

        if strcmpi(order, '1+2') || strcmpi(order, '1')
            ks = round(x(1));
            kh = round(x(2));

            para_xi.kernel.vpara = x(2 + 1: ...
                2 + para_xi.kernel.para_info.vecLength);
            para_alpha.kernel.vpara = x(2 + para_xi.kernel.para_info.vecLength + 1: ...
                2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength);

            if ~isempty(para_xi.kernel.fix)
                para_xi.kernel.vpara(~isnan(para_xi.kernel.vfix)) = para_xi.kernel.vfix(~isnan(para_xi.kernel.vfix));
            end

            if ~isempty(para_alpha.kernel.fix)
                para_alpha.kernel.vpara(~isnan(para_alpha.kernel.vfix)) = para_alpha.kernel.vfix(~isnan(para_alpha.kernel.vfix));
            end

            if return_stat
                para_xi.kernel.vt = model.t(2 + 1: ...
                    2 + para_xi.kernel.para_info.vecLength);
                para_alpha.kernel.vt = model.t(2 + para_xi.kernel.para_info.vecLength + 1: ...
                    2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength);

            end

        end

        if strcmpi(order, '1+2')
            ks = para_alpha.kernel.info.ks;
            kh = para_alpha.kernel.info.kh;
            para_bs.kernel.vpara = x(2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + 1: ...
                2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + para_bs.kernel.para_info.vecLength);

            if ~isempty(para_bs.kernel.fix)
                para_bs.kernel.vpara(~isnan(para_bs.kernel.vfix)) = para_bs.kernel.vfix(~isnan(para_bs.kernel.vfix));
            end

            if return_stat
                para_bs.kernel.vt = model.t(2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + 1: ...
                    2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + para_bs.kernel.para_info.vecLength);
            end

        end

        if strcmpi(order, '2')
            ks = para_alpha.kernel.info.ks;
            kh = para_alpha.kernel.info.kh;
            para_bs.kernel.vpara = x(1: ...
                para_bs.kernel.para_info.vecLength);

            if ~isempty(para_bs.kernel.fix)
                para_bs.kernel.vpara(~isnan(para_bs.kernel.vfix)) = para_bs.kernel.vfix(~isnan(para_bs.kernel.vfix));
            end

            if return_stat
                para_bs.kernel.vt = model.t(1: ...
                    para_bs.kernel.para_info.vecLength);
            end

        end

        if strcmpi(order, '1+2') || strcmpi(order, '2')

            if ~para_bs.use_bssigma
                para_bs.kernel.para.sigma(1, 1) = [para_xi.kernel.para.sigma];
                para_bs.kernel.para.sigma(2, 1) = [para_alpha.kernel.para.sigma];
                para_bs.kernel.lb.sigma = para_bs.kernel.para.sigma;
                para_bs.kernel.ub.sigma = para_bs.kernel.para.sigma;
            end

            if ~para_bs.use_bsnu
                para_bs.kernel.para.nu(1) = [para_xi.kernel.para.nu]; %* 2
                para_bs.kernel.para.nu(2) = [para_alpha.kernel.para.nu]; %* 2
                para_bs.kernel.lb.nu = para_bs.kernel.para.nu;
                para_bs.kernel.ub.nu = para_bs.kernel.para.nu;
            end

            if ~para_xi.use_xid
                para_xi.kernel.para.d = para_alpha.kernel.para.d;
                para_xi.kernel.lb.d = para_alpha.kernel.para.d;
                para_xi.kernel.ub.d = para_alpha.kernel.para.d;
            end

            if ~para_bs.use_bsd
                para_bs.kernel.para.d(1) = [para_xi.kernel.para.d];
                para_bs.kernel.para.d(2) = [para_alpha.kernel.para.d];
                para_bs.kernel.lb.d = para_bs.kernel.para.d;
                para_bs.kernel.ub.d = para_bs.kernel.para.d;
            end

            if ~para_bs.use_bsb
                para_bs.kernel.para.b(1) = [para_xi.kernel.para.b];
                para_bs.kernel.para.b(2) = [para_alpha.kernel.para.b];
                para_bs.kernel.lb.b = para_bs.kernel.para.b;
                para_bs.kernel.ub.b = para_bs.kernel.para.b;
            end

            para_bs.kernel.lb.mu(2) = para_alpha.kernel.para.mu;
            para_bs.kernel.ub.mu(2) = para_alpha.kernel.para.mu;
            para_bs.kernel.para.mu(2) = para_alpha.kernel.para.mu;
        end
    elseif strcmpi(peak_relation, 'free')
        x = x(:);

        if nargin > 5 && ~isempty(model)
            return_stat = true;
        else
            return_stat = false;
        end

        if strcmpi(order, '1+2') || strcmpi(order, '1')
            ks = round(x(1));
            kh = round(x(2));

            para_xi.kernel.vpara = x(2 + 1: ...
                2 + para_xi.kernel.para_info.vecLength);
            para_alpha.kernel.vpara = x(2 + para_xi.kernel.para_info.vecLength + 1: ...
                2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength);

            if ~isempty(para_xi.kernel.fix)
                para_xi.kernel.vpara(~isnan(para_xi.kernel.vfix)) = para_xi.kernel.vfix(~isnan(para_xi.kernel.vfix));
            end

            if ~isempty(para_alpha.kernel.fix)
                para_alpha.kernel.vpara(~isnan(para_alpha.kernel.vfix)) = para_alpha.kernel.vfix(~isnan(para_alpha.kernel.vfix));
            end

            if return_stat
                para_xi.kernel.vt = model.t(2 + 1: ...
                    2 + para_xi.kernel.para_info.vecLength);
                para_alpha.kernel.vt = model.t(2 + para_xi.kernel.para_info.vecLength + 1: ...
                    2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength);

            end

        end

        if strcmpi(order, '1+2')
            ks = para_alpha.kernel.info.ks;
            kh = para_alpha.kernel.info.kh;
            para_bs.kernel.vpara = x(2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + 1: ...
                2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + para_bs.kernel.para_info.vecLength);

            if ~isempty(para_bs.kernel.fix)
                para_bs.kernel.vpara(~isnan(para_bs.kernel.vfix)) = para_bs.kernel.vfix(~isnan(para_bs.kernel.vfix));
            end

            if return_stat
                para_bs.kernel.vt = model.t(2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + 1: ...
                    2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + para_bs.kernel.para_info.vecLength);
            end

        end

        if strcmpi(order, '2')
            ks = para_alpha.kernel.info.ks;
            kh = para_alpha.kernel.info.kh;
            para_bs.kernel.vpara = x(1: ...
                para_bs.kernel.para_info.vecLength);

            if ~isempty(para_bs.kernel.fix)
                para_bs.kernel.vpara(~isnan(para_bs.kernel.vfix)) = para_bs.kernel.vfix(~isnan(para_bs.kernel.vfix));
            end

            if return_stat
                para_bs.kernel.vt = model.t(1: ...
                    para_bs.kernel.para_info.vecLength);
            end

        end

        if strcmpi(order, '1+2') || strcmpi(order, '2')

            kmax_alpha = numel(para_alpha.kernel.para.h);
            alpha_sigma = para_alpha.kernel.para.sigma;
            if isscalar(alpha_sigma), alpha_sigma = repmat(alpha_sigma, kmax_alpha, 1); end
            alpha_nu = para_alpha.kernel.para.nu;
            if isscalar(alpha_nu), alpha_nu = repmat(alpha_nu, kmax_alpha, 1); end
            alpha_d = para_alpha.kernel.para.d;
            if isscalar(alpha_d), alpha_d = repmat(alpha_d, kmax_alpha, 1); end
            alpha_b = para_alpha.kernel.para.b;
            if isscalar(alpha_b), alpha_b = repmat(alpha_b, kmax_alpha, 1); end

            if ~para_bs.use_bssigma
                sigma1 = [para_xi.kernel.para.sigma; alpha_sigma];
                para_bs.kernel.para.sigma = get_ndgrid([sigma1, sigma1], 'list');

                para_bs.kernel.lb.sigma = para_bs.kernel.para.sigma;
                para_bs.kernel.ub.sigma = para_bs.kernel.para.sigma;
            else
                idx_xi = [true(size(para_xi.kernel.para.sigma)); false(size(alpha_sigma))];
                idx_xi = get_ndgrid([idx_xi, idx_xi], 'list');
                sigma1 = [para_xi.kernel.para.sigma; alpha_sigma];
                sigma = get_ndgrid([sigma1, sigma1], 'list');
                para_bs.kernel.para.sigma(idx_xi) = sigma(idx_xi);

                para_bs.kernel.lb.sigma = para_bs.kernel.para.sigma;
                para_bs.kernel.ub.sigma = para_bs.kernel.para.sigma;
            end

            if ~para_bs.use_bsnu
                nu1 = [para_xi.kernel.para.nu; alpha_nu];
                para_bs.kernel.para.nu = get_ndgrid([nu1, nu1], 'list');

                para_bs.kernel.lb.nu = para_bs.kernel.para.nu;
                para_bs.kernel.ub.nu = para_bs.kernel.para.nu;
            else
                idx_xi = [true(size(para_xi.kernel.para.nu)); false(size(alpha_nu))];
                idx_xi = get_ndgrid([idx_xi, idx_xi], 'list');
                nu1 = [para_xi.kernel.para.nu; alpha_nu];
                nu = get_ndgrid([nu1, nu1], 'list');
                para_bs.kernel.para.nu(idx_xi) = nu(idx_xi);
            end

            if ~para_xi.use_xid
                para_xi.kernel.para.d = para_alpha.kernel.para.d;
                para_xi.kernel.lb.d = para_alpha.kernel.para.d;
                para_xi.kernel.ub.d = para_alpha.kernel.para.d;
            end

            if ~para_bs.use_bsd
                d1 = [para_xi.kernel.para.d; alpha_d];
                para_bs.kernel.para.d = get_ndgrid([d1, d1], 'list');

                para_bs.kernel.lb.d = para_bs.kernel.para.d;
                para_bs.kernel.ub.d = para_bs.kernel.para.d;
            else
                idx_xi = [true(size(para_xi.kernel.para.d)); false(size(alpha_d))];
                idx_xi = get_ndgrid([idx_xi, idx_xi], 'list');
                d1 = [para_xi.kernel.para.d; alpha_d];
                d = get_ndgrid([d1, d1], 'list');
                para_bs.kernel.para.d(idx_xi) = d(idx_xi);
            end

            if ~para_bs.use_bsb
                b1 = [para_xi.kernel.para.b; alpha_b];
                para_bs.kernel.para.b = get_ndgrid([b1, b1], 'list');

                para_bs.kernel.lb.b = para_bs.kernel.para.b;
                para_bs.kernel.ub.b = para_bs.kernel.para.b;
            end
            mu1 = [para_xi.kernel.para.mu; para_alpha.kernel.para.mu];
            para_bs.kernel.para.mu = get_ndgrid([mu1, mu1], 'list');
            para_bs.kernel.lb.mu = para_bs.kernel.para.mu;
            para_bs.kernel.ub.mu = para_bs.kernel.para.mu;

        end
    elseif strcmpi(peak_relation, 'harmonic')
        x = x(:);

        if nargin > 5 && ~isempty(model)
            return_stat = true;
        else
            return_stat = false;
        end

        if strcmpi(order, '1+2') || strcmpi(order, '1')
            ks = round(x(1));
            kh = round(x(2));

            para_xi.kernel.vpara = x(2 + 1: ...
                2 + para_xi.kernel.para_info.vecLength);
            para_alpha.kernel.vpara = x(2 + para_xi.kernel.para_info.vecLength + 1: ...
                2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength);

            if ~isempty(para_xi.kernel.fix)
                para_xi.kernel.vpara(~isnan(para_xi.kernel.vfix)) = para_xi.kernel.vfix(~isnan(para_xi.kernel.vfix));
            end

            if ~isempty(para_alpha.kernel.fix)
                para_alpha.kernel.vpara(~isnan(para_alpha.kernel.vfix)) = para_alpha.kernel.vfix(~isnan(para_alpha.kernel.vfix));
            end
            para_alpha.kernel.para.mu = get_harmonic_mu_sigma(para_alpha.kernel.para.mu(ks + 2), 0, ks, kh, para_alpha.no_zero_peak).';

            if return_stat
                para_xi.kernel.vt = model.t(2 + 1: ...
                    2 + para_xi.kernel.para_info.vecLength);
                para_alpha.kernel.vt = model.t(2 + para_xi.kernel.para_info.vecLength + 1: ...
                    2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength);

            end

        end

        if strcmpi(order, '1+2')
            ks = para_alpha.kernel.info.ks;
            kh = para_alpha.kernel.info.kh;
            para_bs.kernel.vpara = x(2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + 1: ...
                2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + para_bs.kernel.para_info.vecLength);

            if ~isempty(para_bs.kernel.fix)
                para_bs.kernel.vpara(~isnan(para_bs.kernel.vfix)) = para_bs.kernel.vfix(~isnan(para_bs.kernel.vfix));
            end

            if return_stat
                para_bs.kernel.vt = model.t(2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + 1: ...
                    2 + para_xi.kernel.para_info.vecLength + para_alpha.kernel.para_info.vecLength + para_bs.kernel.para_info.vecLength);
            end

        end

        if strcmpi(order, '2')
            ks = para_alpha.kernel.info.ks;
            kh = para_alpha.kernel.info.kh;
            para_bs.kernel.vpara = x(1: ...
                para_bs.kernel.para_info.vecLength);

            if ~isempty(para_bs.kernel.fix)
                para_bs.kernel.vpara(~isnan(para_bs.kernel.vfix)) = para_bs.kernel.vfix(~isnan(para_bs.kernel.vfix));
            end

            if return_stat
                para_bs.kernel.vt = model.t(1: ...
                    para_bs.kernel.para_info.vecLength);
            end

        end

        if strcmpi(order, '1+2') || strcmpi(order, '2')

            kmax_alpha = numel(para_alpha.kernel.para.h);
            alpha_sigma = para_alpha.kernel.para.sigma;
            if isscalar(alpha_sigma), alpha_sigma = repmat(alpha_sigma, kmax_alpha, 1); end
            alpha_nu = para_alpha.kernel.para.nu;
            if isscalar(alpha_nu), alpha_nu = repmat(alpha_nu, kmax_alpha, 1); end
            alpha_d = para_alpha.kernel.para.d;
            if isscalar(alpha_d), alpha_d = repmat(alpha_d, kmax_alpha, 1); end
            alpha_b = para_alpha.kernel.para.b;
            if isscalar(alpha_b), alpha_b = repmat(alpha_b, kmax_alpha, 1); end

            if ~para_bs.use_bssigma
                sigma1 = [para_xi.kernel.para.sigma; alpha_sigma];
                para_bs.kernel.para.sigma = get_ndgrid([sigma1, sigma1], 'list');

                para_bs.kernel.lb.sigma = para_bs.kernel.para.sigma;
                para_bs.kernel.ub.sigma = para_bs.kernel.para.sigma;
            else
                idx_xi = [true(size(para_xi.kernel.para.sigma)); false(size(alpha_sigma))];
                idx_xi = get_ndgrid([idx_xi, idx_xi], 'list');
                sigma1 = [para_xi.kernel.para.sigma; alpha_sigma];
                sigma = get_ndgrid([sigma1, sigma1], 'list');
                para_bs.kernel.para.sigma(idx_xi) = sigma(idx_xi);

                para_bs.kernel.lb.sigma = para_bs.kernel.para.sigma;
                para_bs.kernel.ub.sigma = para_bs.kernel.para.sigma;
            end

            if ~para_bs.use_bsnu
                nu1 = [para_xi.kernel.para.nu; alpha_nu];
                para_bs.kernel.para.nu = get_ndgrid([nu1, nu1], 'list');

                para_bs.kernel.lb.nu = para_bs.kernel.para.nu;
                para_bs.kernel.ub.nu = para_bs.kernel.para.nu;
            else
                idx_xi = [true(size(para_xi.kernel.para.nu)); false(size(alpha_nu))];
                idx_xi = get_ndgrid([idx_xi, idx_xi], 'list');
                nu1 = [para_xi.kernel.para.nu; alpha_nu];
                nu = get_ndgrid([nu1, nu1], 'list');
                para_bs.kernel.para.nu(idx_xi) = nu(idx_xi);
            end

            if ~para_xi.use_xid
                para_xi.kernel.para.d = para_alpha.kernel.para.d;
                para_xi.kernel.lb.d = para_alpha.kernel.para.d;
                para_xi.kernel.ub.d = para_alpha.kernel.para.d;
            end

            if ~para_bs.use_bsd
                d1 = [para_xi.kernel.para.d; alpha_d];
                para_bs.kernel.para.d = get_ndgrid([d1, d1], 'list');

                para_bs.kernel.lb.d = para_bs.kernel.para.d;
                para_bs.kernel.ub.d = para_bs.kernel.para.d;
            else
                idx_xi = [true(size(para_xi.kernel.para.d)); false(size(alpha_d))];
                idx_xi = get_ndgrid([idx_xi, idx_xi], 'list');
                d1 = [para_xi.kernel.para.d; alpha_d];
                d = get_ndgrid([d1, d1], 'list');
                para_bs.kernel.para.d(idx_xi) = d(idx_xi);
            end

            if ~para_bs.use_bsb
                b1 = [para_xi.kernel.para.b; alpha_b];
                para_bs.kernel.para.b = get_ndgrid([b1, b1], 'list');

                para_bs.kernel.lb.b = para_bs.kernel.para.b;
                para_bs.kernel.ub.b = para_bs.kernel.para.b;
            end
            mu1 = [para_xi.kernel.para.mu; para_alpha.kernel.para.mu];
            para_bs.kernel.para.mu = get_ndgrid([mu1, mu1], 'list');
            para_bs.kernel.lb.mu = para_bs.kernel.para.mu;
            para_bs.kernel.ub.mu = para_bs.kernel.para.mu;

        end
    end
end