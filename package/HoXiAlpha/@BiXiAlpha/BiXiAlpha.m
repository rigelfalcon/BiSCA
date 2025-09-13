classdef BiXiAlpha < handle & matlab.mixin.Copyable
    %%
    % Authors:
    % Ying Wang
    % Min Li
    % Pedro A. Valdes Sosa
    properties
        f = []
        hos = []
        taperinfo = []
        s = []
        bs = []
        svar
        bsvar
        shatvar
        k0 = []
        kh = []
        ks = []
        khmax = []
        ksmax = []
        kpos = []
        ktype = 'max';

        % xi process
        para_xi = struct();
        % alpha process
        para_alpha = struct();
        % bs para
        para_bs = struct();

        para_fit = struct();

        % para for model compare
        hos_components = 'xi+alpha'; %'xi','alpha' % do we need split xi to two part?
        init_with_subharmonic = 'false'; %'flexible' true false
        region_bs = 'all';

        fx = [];
        fy = [];
        fxy = [];
        fxfy = [];
        nfft = []; % frequency points
        Nf = []; % frequency points used here
        Fs = [];
        T = [];
        df = [];
        width_candidate = 2;

        compact = 'tril'; % 'wedge' 'tril' 'quad1' 'full' %
        multiple_seg = true;
        Ns

        optimization = 'NonlinearOptimization';
        method_xi = 'robust';
        method_alpha = 'maxbs'; %'maxsumbcharmonic'; %'maxbc' %'maxharmonicbc'; %'findmax' 'findpeaks' 'thd' 'maxbcpeak'

        normalization = 'haubrich'; %haubrich skewness hag
        fignormalization = 'haubrich';
        kernel_type = 'harmonic_t_flattop'; %harmonic_t harmonic_t_tapers harmonic_t_flattop
        peak_relation = 'harmonic'; %auto harmonic
        loss = [];

        temp
        % fitting result
        model
        stat
        jac_bc2x

        loss_type = 'bicoherence'; %mse Leonenko BR
        loss_vector = true;
        lh_type = 'bicoherence';
        name = [];
        title = [];
        regularization = 'L1-subharmonic'; %'L1'; %,L1 ,L2[]
        lambda = [0 * 1e-3]; %[0,1e-6]; % [0,1e-3] [1e-6,1e-6]; %[0.001,0.001]; %[0.001,0.001]; %[0,0]; %[0.001,0.01]; %0.01; %0.001

        %Leonenko or BR

        rb = 10; %0 0.5 0.9 1 10
        order = [];
        tfused = [];
        idxused = [];
        tfused_segs = [];
        idxused_segs = [];
        tfused_segs_sep_z = [];
        idxused_segs_sep_z = [];
        iiter = [];
        verbose = false;
        info
        mufix = [];
    end

    properties (Dependent, Hidden)
        x
        lb
        ub

        k
        kmax
        kh2max
        ks2max
        comp
        shat
        smax
        bshat
        bsmax
        sxy
        sxyhat
        sxsysxy
        sxsysxyhat
        bsdenom
        bshatdenom
        bshatvar
        bc
        bchat
        bcmax

        lh
        llh
        nllh
        nllh_s
        nllh_bs
        rss
        mse
        mse_bc
        r2
        r2w
        r2_s
        r2w_s
        r2_bs
        r2w_bs
        r2_bc
        r2w_bc
        r2_bc_th

        dof
        kxi
        N1
        N2
        N1MAX
        N2MAX
        Np
        Np_s
        Np_bs
        Ne

        aic_mse
        bic_mse
        aic_lh
        aic_lh_s
        aic_lh_bs
        bic_lh
        bic_lh_s
        bic_lh_bs
        dht
        dht_n
        summary
    end

    properties (Access = private)
        trend = []
        dn_reg = 0; %1e-4; %1e-4
        s_th = 0;
        p = 1; %0.99; %0.8; % 0.6 0.99;
        q = 1; %0.01; %0.2; % 0.4 0.01;
    end

    methods

        function self = BiXiAlpha(f, s, bs, ks, kh, Fs)

            if nargin > 0

                if nargin < 4 || isempty(ks)
                    ks = 2; %kh has to larger than 2
                end

                if nargin < 5 || isempty(kh)
                    kh = 2; %kh has to larger than 2
                end

                % check and convert all input to double
                f = convert_to_double(f);
                s = convert_to_double(s);
                bs = convert_to_double(bs);
                ks = convert_to_double(ks);
                kh = convert_to_double(kh);
                Fs = convert_to_double(Fs);

                if nargin < 5 || isempty(Fs)
                    Fs = 2 * max(f);
                end

                self.f = f;
                self.s = s;
                self.bs = bs;
                self.ks = ks;
                self.kh = kh;
                self.Fs = Fs;
                self.ksmax = ks;
                self.khmax = kh;
            end

        end

        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function set.s(self, val)

            if any(val < 0)
                error('spectrum must be positive')
            end

            self.s = val;
        end

        function set.hos(self, val)
            if isobject(val)
                warning('off', 'MATLAB:structOnObject');
                self.hos = struct(val);
                warning('on', 'MATLAB:structOnObject');
            else
                self.hos = val;
            end

            if ~isempty(val)

                fields_remove = {'W2', 'W3', 'bs', 'bc', 'fy', 'fx', 'X'};
                fields = fieldnames(self.hos);
                idx_loc = endsWith(fields, 'loc');
                fields_remove = [fields_remove, fields(idx_loc)'];
                self.hos = set_field(self.hos, fields_remove, []);

                self.normalization = self.hos.normalization;

                self.svar = self.hos.svar; %#ok<*MCSUP> % self.shatvar = self.hos.svar;
                self.hos.svar = [];
                self.shatvar = [];

                self.bsvar = self.hos.bsvar;
                self.hos.bsvar = [];
                self.bsvar = tril2full(self.bsvar, 0, false);
                self.bsvar(self.bsvar == 0) = NaN;

                self.taperinfo = self.hos;
                self.taperinfo.bsdenom = []; % save memory
                self.taperinfo = update_taperinfo(self.taperinfo);
                self.Ns = size(self.s, 2);

                if self.Ns ~= size(self.bs, 3)
                    error('the number of segment is not match')
                end

            end

        end

        function set.f(self, val)
            self.f = val;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% get
        function get_bifrequency(self)
            [self.fx, self.fy] = ndgrid(self.f);
            self.fxy = reshape(self.fx + self.fy, self.Nf * self.Nf, 1);
            self.para_bs.idx_fxfy = logical(mat2tril(mat2atriu(true(size(self.fx)), 0, false), 0, false));
            self.fxfy = [self.fx(self.para_bs.idx_fxfy), self.fy(self.para_bs.idx_fxfy)];
        end

        function ks2max = get.ks2max(self)
            ks2max = self.ksmax - self.ks;
        end

        function kmax = get.kmax(self)
            kmax = self.k0 + self.ksmax + self.khmax;
        end

        function kmax = get.k(self)
            kmax = self.k0 + self.ks + self.kh;
        end

        function T = get.T(self)
            T = self.nfft / self.Fs;
        end

        function Nf = get.Nf(self)
            Nf = length(self.f);
        end

        function nfft = get.nfft(self)
            nfft = round(self.Fs / self.df);
        end

        function df = get.df(self)
            fmax = max(self.f);
            df = fmax / length(self.f);
        end

        function shat = get.shat(self)
            shat = pred_s(self.f, self.ks, self.kh, self.para_xi, self.para_alpha);
        end

        function bshat = get.bshat(self)
            bshat = pred_bs(self.f, self.para_bs, self.fxfy);
        end

        function sxy = get.sxy(self)

            if ~isempty(self.bs)
                opt.ExtrapolationMethod = 'none';
                sxy = get_sxy(self.s, self.f, opt);
                sxy(isnan(sxy)) = 0;
            end

        end

        function sxyhat = get.sxyhat(self)

            if ~isempty(self.bs)
                opt.ExtrapolationMethod = 'none';
                sxyhat = get_sxy(self.shat, self.f, opt);
                sxyhat(isnan(sxyhat)) = 0;
            end

        end

        function sxsysxy = get.sxsysxy(self)

            if ~isempty(self.bs)

                if ~isempty(self.fxy)
                    sxsysxy = get_sxsysxy(self.f, self.s, [], self.fxy);
                else
                    sxsysxy = get_sxsysxy(self.f, self.s);
                end

            end

        end

        function sxsysxyhat = get.sxsysxyhat(self)

            if ~isempty(self.bs)

                if ~isempty(self.fxy)
                    sxsysxyhat = get_sxsysxy(self.f, self.shat, [], self.fxy);
                else
                    sxsysxyhat = get_sxsysxy(self.f, self.shat);
                end

            end

        end

        function bsdenom = get.bsdenom(self)

            if ~isempty(self.bs)

                switch self.normalization
                    case {'hagihira', 'hag'}
                        bsdenom = self.hos.bsdenom;
                    otherwise
                        bsdenom = sqrt(self.sxsysxy);

                end

            end

        end

        function bshatdenom = get.bshatdenom(self)

            if ~isempty(self.bs)

                switch self.normalization
                    case {'hagihira', 'hag'}
                        bshatdenom = self.hos.bsdenom;
                    otherwise

                        if isfield(self.para_bs, 'bshatdenom_pre') && ~isempty(self.para_bs.bshatdenom_pre)
                            bshatdenom = self.para_bs.bshatdenom_pre;
                        else
                            bshatdenom = sqrt(self.sxsysxyhat);
                        end

                end

            end

        end

        function bshatvar = get.bshatvar(self)

            if ~isempty(self.bs)

                bshatvar = calc_bshatvar(self, self.bshatdenom .^ 2);
            end

        end

        function bshatvar = calc_bshatvar(self, bshatvar_unscale)

            if ~isempty(self.bs)

                if ~isempty(self.bsvar)

                    if nargin < 2
                        bshatvar_unscale = self.sxsysxyhat;
                    end

                    bshatvar_unscale(bshatvar_unscale == 0) = NaN;
                    ratio = mean(self.bsvar ./ (bshatvar_unscale), [1, 2], 'omitnan'); %sqrt
                    bshatvar = bshatvar_unscale .* ratio;
                else
                    bshatvar = bshatvar_unscale;
                end

            end

        end

        function bc = get.bc(self)
            bc = calc_bc(self, self.bs, self.bsdenom);
        end

        function bchat = get.bchat(self)
            bchat = calc_bc(self, self.bshat, self.bshatdenom);
        end

        function x = get.x(self)

            switch self.order
                case {'1'}
                    x = [self.vec_para('k'); self.vec_para('xi'); self.vec_para('alpha')];
                case {'2'}
                    x = [self.vec_para('bs')];
                case {'1+2'}
                    x = [self.vec_para('k'); self.vec_para('xi'); self.vec_para('alpha'); self.vec_para('bs')];
                otherwise
                    error('no such order')
            end

        end

        function lb = get.lb(self)

            switch self.order
                case '1'
                    lb = [bound_k(self); bound_xi(self); bound_alpha(self)];
                case '2'
                    lb = [bound_bs(self)];
                case '1+2'
                    lb = [bound_k(self); bound_xi(self); bound_alpha(self); bound_bs(self)];
            end

        end

        function ub = get.ub(self)

            switch self.order
                case '1'
                    [~, ubk] = bound_k(self);
                    [~, ubxi] = bound_xi(self);
                    [~, ubalpha] = bound_alpha(self);
                    ub = [ubk; ubxi; ubalpha];
                case '2'
                    [~, ub] = bound_bs(self);
                case '1+2'
                    [~, ubk] = bound_k(self);
                    [~, ubxi] = bound_xi(self);
                    [~, ubalpha] = bound_alpha(self);
                    [~, ubbs] = bound_bs(self);
                    ub = [ubk; ubxi; ubalpha; ubbs];
            end

        end

        function kxi = get.kxi(self)

            if strcmpi(self.para_bs.type, 'bigaussian+xi+delta') || strcmpi(self.para_bs.type, 'trigaussian+xi+delta')
                kxi = 2;
            elseif strcmpi(self.para_bs.type, 'bigaussian+xi') || strcmpi(self.para_bs.type, 'trigaussian+xi') || strcmpi(self.para_bs.type, 'tstudent')
                kxi = 1;
            else
                error('no such method')
            end

        end

        function N1 = get.N1(self)

            N1 = self.k + self.kxi;
        end

        function N2 = get.N2(self)
            N2 = self.N1 ^ 2;
        end

        function N1MAX = get.N1MAX(self)

            N1MAX = self.kmax + self.kxi;
        end

        function N2MAX = get.N2MAX(self)
            N2MAX = self.N1MAX ^ 2;
        end

        function Np = get.Np(self)
            self.order = '1+2';
            x0 = self.x;
            ub = self.ub;
            lb = self.lb;
            x0 = remove_fixed_para(x0, lb, ub);
            Np = length(x0); %put zero as none active parameter

        end

        function Np_s = get.Np_s(self)
            self.order = '1';
            x0 = self.x;
            ub = self.ub;
            lb = self.lb;
            x0 = remove_fixed_para(x0, lb, ub);
            Np_s = length(x0); %put zero as none active parameter
        end

        function Np_bs = get.Np_bs(self)
            self.order = '2';
            x0 = self.x;
            ub = self.ub;
            lb = self.lb;
            x0 = remove_fixed_para(x0, lb, ub);
            Np_bs = length(x0); %put zero as none active parameter
        end

        function [t, ci, se] = calc_bc_t(self, component, ispolar)
            jac_r2x = self.model.jacobian;

            jac_bc2x = calc_jac_bc2x(self, component, ispolar);

            bchat = self.bchat(self.tfused);

            if ispolar
                bchat = [abs(bchat); angle(bchat)];
            else
                bchat = [real(bchat); imag(bchat)];
            end

            [ci, se] = nlpar_ci_se_bc(bchat, self.model.residual, jac_r2x, jac_bc2x);
            t = bchat ./ se;
            Nbc = length(bchat) / 2;

            % t value of angle nolonger has the meaning of phase, therefore, here we use complex form
            t = t(1:Nbc) + 1i * t(Nbc + 1:2 * Nbc);
            ci = ci(1:Nbc) + 1i * ci(Nbc + 1:2 * Nbc);
            se = se(1:Nbc) + 1i * se(Nbc + 1:2 * Nbc);
            t = bmap(t, self.tfused, self.bchat);
            ci = bmap(ci, self.tfused, self.bchat);
            se = bmap(se, self.tfused, self.bchat);

        end

        function jac_bc2x = calc_jac_bc2x(self, component, ispolar)
            % since the sample is not fixed for each bootstrap, we need to estimate the jacobian
            tic
            fun_x2bc = @(x) self.x2bc(component, x, ispolar);
            para = self.x(~self.para_fit.parainfo.xIndices.fixed);
            [jac_bc2x, err] = jacobianest(fun_x2bc, para);
            jac_bc2x = sparse(jac_bc2x);
            t = toc;
            disp(['[calc_jac_bc2x] time cost for estimation = ', num2str(t), ' s']);
        end

        function [bchat] = x2bc(self, component, para, ispolar)
            % only active parameters
            para = recover_fixed_para(para, [], [], self.para_fit.parainfo);
            [ks, kh, para_xi_perturb, para_alpha_perturb, para_bs_perturb] = vecpara2para(para, self.order, self.para_xi, self.para_alpha, self.para_bs);
            %% this is for the bootstrap, so only consider the linear parameter height
            para_xi = self.para_xi;

            para_alpha = self.para_alpha;

            para_bs = self.para_bs;
            para_bs.kernel.para.h = para_bs_perturb.kernel.para.h;

            bshat = bshat_component(self, component, para_bs);

            switch self.normalization
                case {'hagihira', 'hag'}
                    bshatdenom = self.hos.bsdenom; %#ok<*PROPLC>
                otherwise
                    shat = pred_s(self.f, ks, kh, para_xi, para_alpha);
                    sxsysxyhat = get_sxsysxy(self.f, shat, [], self.fxy);
                    bshatdenom = sqrt(sxsysxyhat);
            end

            bchat = calc_bc(self, bshat, bshatdenom);
            bchat = bchat(self.tfused);

            if ispolar
                bchat = [abs(bchat); angle(bchat)];
            else
                bchat = [real(bchat); imag(bchat)];
            end

        end

        function dof = get.dof(self)
            dof = self.Np;
        end

        function mse = get.mse(self)
            bs = self.bs; %#ok<*PROP> %bs(isnan(bs)) = 0;
            bs = compact_bs(bs, self.idxused);

            bshat = self.bshat; %bshat(isnan(bshat)) = 0;
            bshat = compact_bs(bshat, self.idxused);
            mse = mean((self.shat - self.s) .^ 2) + mean((abs(bs - bshat) .^ 2), "all");
        end

        function mse = get.mse_bc(self)
            bc = self.bc; %#ok<*PROP> %bs(isnan(bs)) = 0;
            bc = compact_bs(bc, self.idxused);

            bchat = self.bchat; %bshat(isnan(bshat)) = 0;
            bchat = compact_bs(bchat, self.idxused);
            mse = mean((self.shat - self.s) .^ 2) + mean((abs(bc - bchat) .^ 2), "all");
        end

        function rss = get.rss(self)
            bs = self.bs; %bs(isnan(bs)) = 0;
            bs = compact_bs(bs, self.idxused);

            bshat = self.bshat; %bshat(isnan(bshat)) = 0;
            bshat = compact_bs(bshat, self.idxused);
            rss = sum((self.shat - self.s) .^ 2) + sum((abs(bs - bshat) .^ 2), "all");
        end

        function r2 = get.r2(self)
            Ns = self.Ns;
            bs = reshape(self.bs, [], Ns);
            y = [self.s; real(bs); imag(bs)];
            bshat = reshape(self.bshat, [], 1);
            yhat = [self.shat; real(bshat); imag(bshat)];
            r2 = calc_goodnessOfFit(y, yhat);
        end

        function r2w = get.r2w(self)
            Ns = self.Ns;
            bs = reshape(self.bs, [], Ns);
            y = [log(self.s); real(bs); imag(bs)];
            bshat = reshape(self.bshat, [], 1);
            yhat = [log(self.shat); real(bshat); imag(bshat)];
            bsvar = reshape(self.bshatdenom, [], Ns); %.^ 2
            yvar = [self.svar; bsvar; bsvar];
            r2w = calc_r2w(y, yhat, yvar);
        end

        function r2 = get.r2_s(self)
            y = [self.s];
            yhat = [self.shat];
            r2 = calc_goodnessOfFit(y, yhat);
        end

        function r2w = get.r2w_s(self)
            y = [self.s];
            yhat = [self.shat];
            yvar = [self.svar];
            r2w = calc_r2w(y, yhat, yvar);
        end

        function r2 = get.r2_bs(self)
            Ns = self.Ns;
            bs = reshape(self.bs, [], Ns);
            y = [real(bs); imag(bs)];
            bshat = reshape(self.bshat, [], 1);
            yhat = [real(bshat); imag(bshat)];
            r2 = calc_goodnessOfFit(y, yhat);
        end

        function r2 = get.r2_bc(self)
            Ns = self.Ns;
            bc = reshape(self.bc, [], Ns);
            y = [real(bc); imag(bc)];
            bchat = reshape(self.bchat, [], 1);
            yhat = [real(bchat); imag(bchat)];
            r2 = calc_goodnessOfFit(y, yhat);
        end

        function r2 = get.r2_bc_th(self)
            Ns = self.Ns;
            bc = reshape(self.bc(~isnan(self.stat.raw.max.bcth)), [], Ns);
            y = [real(bc); imag(bc)];
            bchat = reshape(self.bchat(~isnan(self.stat.raw.max.bcth)), [], 1);
            yhat = [real(bchat); imag(bchat)];
            r2 = calc_goodnessOfFit(y, yhat);
        end

        function r2w = get.r2w_bs(self)
            Ns = self.Ns;
            bs = reshape(self.bs, [], Ns);
            y = [real(bs); imag(bs)];
            bshat = reshape(self.bshat, [], 1);
            yhat = [real(bshat); imag(bshat)];
            bsvar = reshape(self.bshatdenom, [], Ns); %.^ 2
            yvar = [bsvar; bsvar];
            r2w = calc_r2w(y, yhat, yvar);
        end

        function r2w = get.r2w_bc(self)
            Ns = self.Ns;
            bc = reshape(self.bc, [], Ns);
            y = [real(bc); imag(bc)];
            bchat = reshape(self.bchat, [], 1);
            yhat = [real(bchat); imag(bchat)];
            bcvar = reshape(self.bchatdenom, [], Ns); %.^ 2
            yvar = [bcvar; bcvar];
            r2w = calc_r2w(y, yhat, yvar);
        end

        function lh = get.lh(self)
            lh = exp(self.llh);
        end

        function Ne = get.Ne(self)
            % number of elements of data
            Ne = self.Nf + (self.Nf * self.Nf + self.Nf) / 2;
        end

        function llh = get.llh(self)
            % the larger llh is the better model is
            llh = -self.nllh;
        end

        % Ns:Ns*NW*Nf
        function aic_mse = get.aic_mse(self)
            % the larger aic_mse is the worse model is
            aic_mse = 2 * self.Np + self.Ne * log(self.mse);
        end

        function bic_mse = get.bic_mse(self)
            n = self.Nf * self.Ns * round(self.hos.NW * 2 - 1);
            bic_mse = self.Np * log(n) + n * log(self.mse);
        end

        function aic_lh_s = get.aic_lh_s(self)
            % warning('seperate complex')
            aic_lh_s = 2 * self.Np_s + 2 * self.nllh_s;
        end

        function aic_lh_bs = get.aic_lh_bs(self)
            % warning('seperate complex')
            aic_lh_bs = self.Np_bs + 2 * self.nllh_bs;
        end

        function aic_lh = get.aic_lh(self)
            % the larger aic_lh is the worse model is
            aic_lh = 2 * self.Np + 2 * self.nllh;
        end

        function bic_lh = get.bic_lh(self)
            % the larger bic is the worse model is
            n = self.Nf * self.Ns * round(self.hos.NW * 2 - 1);
            bic_lh = self.Np * log(n) + 2 * self.nllh;
        end

        function bic_lh_s = get.bic_lh_s(self)
            % the larger bic is the worse model is

            n = self.Nf * self.Ns * round(self.hos.NW * 2 - 1);
            bic_lh_s = self.Np_s * log(n) + 2 * self.nllh_s;
        end

        function bic_lh_bs = get.bic_lh_bs(self)
            % the larger bic is the worse model is
            n = self.Nf * self.Ns * round(self.hos.NW * 2 - 1);
            bic_lh_bs = self.Np_bs * log(n) + 2 * self.nllh_bs;
        end

        function nllh = get.nllh(self)

            bshatvar = self.bshatdenom;
            bc = self.bc;
            nllh = sum(abs(calc_loss(self, self.lh_type, self.s, self.shat, self.svar, bc, self.bshat, [], bshatvar)) .^ 2);
            if imag(nllh) ~= 0
                nllh = nan;
            end

        end

        function nllh_s = get.nllh_s(self)
            nllh_s = sum(abs(calc_loss_s(self, self.lh_type, self.s, self.shat, self.svar)) .^ 2);

            if imag(nllh_s) ~= 0
                nllh_s = nan;
            end

        end

        function nllh_bs = get.nllh_bs(self)

            switch self.loss_type
                case {'bicoherence'}
                    bsvar = self.bsdenom;
                    bshatvar = self.bshatdenom;
                case {'HeteroEstVar', 'HeteroscedasticityEstVar'}
                    bsvar = self.bsvar;
                    bshatvar = self.bshatvar;
            end

            bc = self.bc;
            nllh_bs = sum(abs(calc_loss_bs(self, self.lh_type, bc, self.bshat, bsvar, bshatvar)) .^ 2);

            if imag(nllh_bs) ~= 0
                nllh_bs = nan;
            end

        end

        function dht = get.dht(self)
            [~, scomp] = pred_s_alpha(self.f, self.kh, self.para_alpha);
            dht = 10 * log10(sum(scomp(:, 2:end), 'all') / sum(scomp(:, 1)));
        end

        function dht_n = get.dht_n(self)
            [~, scomp_alpha] = pred_s_alpha(self.f, self.kh, self.para_alpha);
            [~, scomp_xi] = pred_s_xi(self.f, self.para_xi);
            dht_n = 10 * log10((sum(scomp_xi, 'all') + sum(scomp_alpha(:, 2:end), 'all')) / sum(scomp_alpha(:, 1)));
        end

        function summary = get.summary(self)

            if isempty(self.title)

                if strcmpi(self.hos_components, 'xi+alpha')
                    summary.title = "xi+alpha (unrestricted, fit xi and alpha in bispectrum)";
                elseif strcmpi(self.hos_components, 'alpha')
                    summary.title = "alpha (restricted, drop xi in bispectrum)";
                elseif strcmpi(self.hos_components, 'xi')
                    summary.title = "xi (restricted, drop alpha in bispectrum)";
                else
                    error('no such method')
                end

            end

            summary.ks = self.ks;
            summary.kh = self.kh;
            summary.num_para = self.Np;
            summary.num_para_s = self.Np_s;
            summary.num_para_bs = self.Np_bs;
            summary.loss = self.loss;
            summary.nllh_s = self.nllh_s;
            summary.nllh_bs = self.nllh_bs;
            summary.aic_lh = self.aic_lh;
            summary.aic_lh_s = self.aic_lh_s;
            summary.aic_lh_bs = self.aic_lh_bs;
            summary.bic_lh = self.bic_lh;
            summary.bic_lh_s = self.bic_lh_s;
            summary.bic_lh_bs = self.bic_lh_bs;

            summary.r2 = self.r2;
            summary.mse = self.mse;
            summary.aic_mse = self.aic_mse;
            summary = struct2table(summary);
        end

    end

end