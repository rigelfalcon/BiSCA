function loss = ObjFunBiXiAlpha(self, para, record)

    if self.iiter == 0
        self.temp.para = NaN(length(para), 2500);
    else
        self.temp.para(:, int64(self.iiter + 1)) = para;
    end

    if nargin < 3
        record = [];
    end

    % record: avoid redundant calculation
    % get para to fit from vectorized para
    para = recover_fixed_para(para, [], [], self.para_fit.parainfo);
    [ks, kh, para_xi, para_alpha, para_bs] = vecpara2para(para, self.order, self.para_xi, self.para_alpha, self.para_bs, [], self.peak_relation);

    %% spectrum
    if strcmpi(self.order, '1')
        shat = pred_s(self.f, ks, kh, para_xi, para_alpha);
    elseif strcmpi(self.order, '1+2')
        shat = pred_s(self.f, ks, kh, para_xi, para_alpha);
    end

    %% bispectrum
    if strcmpi(self.order, '1+2')
        bshat = pred_bs(self.f, para_bs, self.fxfy);

        switch self.loss_type
            case {'bicoherence'}

                bshatdenom = sqrt(get_sxsysxy(self.f, shat, [], self.fxy)); %sxsysxyhat
                bshatdenom(bshatdenom == 0) = NaN;

                bc = record.bc;
                bshatvar = bshatdenom;

            otherwise
        end

        if strcmpi(self.region_bs, "mu+diag")

            if self.para_alpha.isubharmonic %&& (strcmpi(self.hos_components,'xi+alpha') || strcmpi(self.hos_components, 'alpha'))
                idx_mu = find_closest(self.f, self.para_alpha.kernel.para.mu * 2);
            else
                idx_mu = find_closest(self.f, self.para_alpha.kernel.para.mu);
            end

            bshat = [bshat(idx_mu, :)'; bshat(record.idx_diag)];
            bshatvar = [bshatvar(idx_mu, :)'; record.bshatvar_diag];
        end

    elseif strcmpi(self.order, '2')
        bshat = pred_bs(self.f, para_bs, self.fxfy);

        if strcmpi(self.region_bs, "boundary")
            bshat = bshat(record.idx_mu, :)';
            bshatvar = record.bshatvar_boundary;
        elseif strcmpi(self.region_bs, "diag")
            bshat = bshat(record.idx_diag);
            bshatvar = record.bshatvar_diag;
        elseif strcmpi(self.region_bs, "offdiag") || strcmpi(self.region_bs, "all") || isempty(self.region_bs)

            bc = record.bc;
            bshatvar = record.bshatvar;
        else
            error('not support')
        end

    end

    %% lost function
    if strcmpi(self.order, '1')
        %%{
        loss = calc_loss_s(self, self.loss_type, self.s, shat, self.svar); %sep_complex
        %%}
    elseif strcmpi(self.order, '2')
        loss = calc_loss_bs(self, self.loss_type, bc, bshat, [], bshatvar); %sep_complex
    elseif strcmpi(self.order, '1+2')
        %%{
        [loss] = calc_loss(self, self.loss_type, self.s, shat, self.svar, bc, bshat, [], bshatvar);
        %%}
    end

    penalty = calc_penalty(self, para_xi, para_alpha, para_bs);
    loss = [loss; penalty];

    self.iiter = self.iiter + 1;

    if strcmpi(self.order, '1')
        itv_show = 50;

        if self.verbose && (mod(self.iiter, itv_show) == 1 || self.iiter == 1)
            s = meanlog10(self.s, 2);
            figure(10);
            plot(self.f, log10([s, shat]))
            drawnow
            sgtitle(['iter=', num2str(self.iiter)])
            if self.iiter == 1
            end

        end

    elseif strcmpi(self.order, '2')
        itv_show = 400;
        if self.verbose && (mod(self.iiter, itv_show) == 1 || self.iiter == 1)
            figure(11)

            bchat = real(bshat ./ mean(self.bshatdenom, 3));
            subplot(2, 1, 1)
            surfbc(self.fx, self.fy, real(mean(self.bc, 3)));
            clim('auto')
            subplot(2, 1, 2)
            surfbc(self.fx, self.fy, bchat);
            clim('auto')
            drawnow
            sgtitle(['iter=', num2str(self.iiter)])
            if self.iiter == 1
            end

        end

    elseif strcmpi(self.order, '1+2')
        itv_show = 400;
        if self.verbose && (mod(self.iiter, itv_show) == 1 || self.iiter == 1)
            figure(12)
            subplot(3, 1, 1)
            plot(log10([self.s, shat]))

            subplot(3, 2, 3)
            bc = mean(self.bc, 3);
            surfbc(self.fx, self.fy, real(bc));
            clim('auto')

            subplot(3, 2, 4)
            surfbc(self.fx, self.fy, imag(bc));
            clim('auto')

            subplot(3, 2, 5)
            bchat = (bshat ./ mean(self.bshatdenom, 3));
            surfbc(self.fx, self.fy, real(bchat));
            clim('auto')
            drawnow
            subplot(3, 2, 6)
            surfbc(self.fx, self.fy, imag(bchat));
            clim('auto')
            drawnow
            sgtitle(['iter=', num2str(self.iiter)])
            if self.iiter == 1
            end

        end

    end

    self.model.loss_history(self.iiter, :) = norm(loss(~isnan(loss)));
    loss = convert_to_double(loss);

end