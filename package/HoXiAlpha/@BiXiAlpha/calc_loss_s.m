function loss = calc_loss_s(self, loss_type, s, shat, svar, shatvar)

    if isempty(loss_type)
        loss_type = self.loos_type;
    end

    if nargin < 6 || isempty(shatvar)
        shatvar = self.shatvar;
    end

    if nargin < 5 || isempty(svar)
        svar = self.svar;
    end

    switch loss_type
        case {'bicoherence'}
            loss = ((log(shat) - log(s)) ./ sqrt(svar));
        case {'Hetero', 'Heteroscedasticity', 'HeteroEstVar', 'HeteroscedasticityEstVar'}
            shat(shat < eps) = eps;
            loss = ((log(shat) - log(s)) ./ sqrt(svar)); % log transform make homosdestic
        case 'Leonenko'
            loss = ((shat - s) ./ shat) .^ 2;
        case 'mse'
            loss = (shat - s) .^ 2;
        case {'BR', 'AIC', 'MCE'} %,'BRabs'
            idx_calc = ~isnan(shat) & shat ~= 0;
            loss = log(shat(idx_calc)) + s(idx_calc) ./ shat(idx_calc); %((2*pi)^(-2)/self.Nf)*
        case {'BRNg', 'BrillingerNongaussian'} %,'BRabs'
            idx_calc = ~isnan(shat) & shat ~= 0;
            loss = ((s(idx_calc) - shat(idx_calc)) .^ 2) ./ shat(idx_calc);
        case {'Rice'}
            idx_calc = ~isnan(shat) & shat ~= 0;
            loss = ((s(idx_calc) - shat(idx_calc)) ./ shat(idx_calc)) .^ 2;
        otherwise
            error('no such loss_type')
    end

    loss = loss(:);
    loss(isnan(loss) | isinf(loss)) = inf;

    if strcmpi(self.para_fit.optfun, 'fmincon')
        loss = sum(loss, "all");
    end

    loss = convert_to_double(loss);

end