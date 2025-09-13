function loss = calc_loss_bs(self, loss_type, bs, bshat, bsvar, bshatvar)

    if isempty(loss_type)
        loss_type = self.loos_type;
    end

    idx_calc = self.tfused_segs;

    switch loss_type
        case {'bicoherence'}
            bc = bs;
            bshatdenom = bshatvar + self.dn_reg; %sqrt
            loss = cat(3, (real(bc) - real(bshat) ./ bshatdenom), ...
                (imag(bc) - imag(bshat) ./ bshatdenom)); %abs
            idx_calc = self.tfused_segs_sep_z;
        case {'bicoherence_complecated'}

            bsdenom = bsvar + self.dn_reg; %sqrt
            bshatdenom = bshatvar + self.dn_reg; %sqrt
            loss = cat(3, (real(bs) ./ bsdenom - real(bshat) ./ bshatdenom), ...
                (imag(bs) ./ bsdenom - imag(bshat) ./ bshatdenom)); %abs
            idx_calc = self.tfused_segs_sep_z;
        case {'HeteroEstVar', 'HeteroscedasticityEstVar'}

            if isempty(bshatvar)
                warning('bshatvar is empty, use bsvar')
                bshatvar = bsvar;
            end

            bshatstd = sqrt(bshatvar) + self.dn_reg;
            loss = cat(3, (real(bs) - real(bshat)) ./ bshatstd, (imag(bs) - imag(bshat)) ./ bshatstd);
            idx_calc = self.tfused_segs_sep_z;

        case {'Hetero', 'Heteroscedasticity'}
            bsstd = sqrt(bsvar) + self.dn_reg;
            loss = cat(3, (real(bs) - real(bshat)) ./ bsstd, (imag(bs) - imag(bshat)) ./ bsstd);
            idx_calc = self.tfused_segs_sep_z;
        case {'Leonenko'}
            loss = (abs(bs(idx_calc) - bshat(idx_calc)) .^ 2) ./ (bshatvar(idx_calc));
        case {'mse'}
            loss = (abs(bshat(idx_calc) - bs(idx_calc)) .^ 2);
        case {'BR', 'AIC', 'BRNg'}

            switch lower(self.normalization)
                case {'haubrich'}
                    loss = (abs(bs(idx_calc) - bshat(idx_calc)) .^ 2) ./ (bshatvar(idx_calc));
                case {'skewness'}
                    loss = (abs(bs(idx_calc) - bshat(idx_calc)) .^ 2) ./ (bshatvar(idx_calc)); %+1e-10 (2*pi/self.Nf)*
                otherwise
                    error('no such method')
            end

        case {'Rice'}
            idx_calc = idx_calc & bshat ~= 0 & bshat > (10 * eps);
            loss = abs((bs(idx_calc) - bshat(idx_calc)) ./ (bshat(idx_calc))) .^ 2; %any(isinf(loss))||any(isnan(loss))
            loss(isinf(loss)) = 0;
        case {'MCE'}
            loss = log(abs(bshat(idx_calc)) .^ 2) + abs(bs(idx_calc)) .^ 2 ./ abs(bshat(idx_calc)) .^ 2;
        otherwise
            error('no such loss_type')
    end

    loss = loss(idx_calc);
    loss(isnan(loss) | (isinf(loss) & loss > 0)) = inf;

    if strcmpi(self.para_fit.optfun, 'fmincon')
        loss = sum(loss, "all");
    end

    loss = convert_to_double(loss);

end