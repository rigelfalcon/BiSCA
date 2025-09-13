classdef BiXiAlphaBootstrap < handle & matlab.mixin.Copyable

    properties
        x = [];
        hos
        multiple_seg = false;
        xialpha = BiXiAlpha();
        nlng
        null_xi
        null_alpha
        Nb = 200; %100 50
        Ns = []
        Ns_eff = []
        khall
        ksall
        metric = 'nllh'; %nllh aic_lh bic_lh
        logging
        alpha_level = 0.001;
        p_threshold = 0.001;
        crt
        idx_sample
        rngseed
    end

    properties (Hidden)
        verbose = []
        tf_runtest = true %true ~ispc
        models
        modelsnull
        shatres
        bshatres
        maxworker = 25;
        maxattemp = 50;
    end

    properties (Dependent, Hidden)
        summary

    end

    methods

        function self = BiXiAlphaBootstrap(hos, ksall, khall)

            if nargin > 0

                if nargin > 1
                    self.hos = hos;
                    self.ksall = ksall;
                    self.khall = khall;
                elseif nargin == 1
                    self.hos = hos;
                    self.ksall = 0:2;
                    self.khall = 1:8;
                end

                self.Ns = self.hos.Ns;
            end

        end

        function init(self)

            if isnumeric(self.hos.name)
                self.hos.name = num2str(self.hos.name);
            end

            self.hos.s; %excute calculation

            % copy hos compact meta info to xialpha
            if isobject(self.hos) || ~self.hos.update_lock
                self.x = self.hos.X;
                self.hos.compact_result;
                warningState = warning('off');
                self.hos = struct(self.hos);
                warning(warningState);
                hos = self.hos; %#ok<*PROP>
                hos.sres = [];
                hos.bsres = [];
            end

            % generate multiple samples for the regression
            if self.multiple_seg

                if self.hos.tf_return_seg
                    s = self.hos.s;
                    bs = self.hos.bs;
                elseif ~self.hos.tf_return_seg && self.hos.tf_return_res
                    s = exp(log(self.hos.s) + self.hos.sres); %#ok<*CPROP>
                    bs = self.hos.bs + self.hos.bsres;
                else
                    error('need return seg or keep res')
                end

                if ~isempty(self.Ns) && self.Ns ~= self.hos.Ns
                    self.idx_sample = randi(self.hos.Ns, [self.Ns, 1]);
                    s = s(:, self.idx_sample);
                    bs = bs(:, :, self.idx_sample);
                    hos.bsdenom = self.hos.bsdenom(:, :, self.idx_sample);
                end

            else
                s = self.hos.s;
                bs = self.hos.bs;
            end

            for iks = length(self.ksall):-1:1
                ks = self.ksall(iks);

                for ikh = length(self.khall):-1:1
                    kh = self.khall(ikh);
                    self.xialpha(ikh, iks) = BiXiAlpha(hos.f, s, bs, ks, kh, hos.Fs);
                    self.xialpha(ikh, iks).hos = hos;
                    self.xialpha(ikh, iks).ksmax = ks;
                    self.xialpha(ikh, iks).khmax = kh;

                    self.xialpha(ikh, iks).hos_components = 'xi+alpha'; % 'xi+alpha' 'alpha' 'xi'

                    self.xialpha(ikh, iks).loss_type = 'bicoherence'; %BR  Leonenko BRNg Rice Hetero HeteroEstVar bicoherence
                    self.xialpha(ikh, iks).lh_type = 'bicoherence'; %BR  Leonenko Hetero HeteroEstVar
                    self.xialpha(ikh, iks).method_alpha = 'findmax';
                    self.xialpha(ikh, iks).peak_relation = 'free';
                    self.xialpha(ikh, iks).kernel_type = 'harmonic_t_free';

                    self.xialpha(ikh, iks).regularization = 'elastic_constraint_mu_sigma_chi_xi';
                    self.xialpha(ikh, iks).para_bs.no_xi_alpha = true;
                    self.xialpha(ikh, iks).para_bs.use_bssigma = false;
                    self.xialpha(ikh, iks).para_bs.use_bsnu = false;
                    self.xialpha(ikh, iks).para_bs.use_bsd = false;
                    self.xialpha(ikh, iks).para_bs.use_bsb = false;
                    self.xialpha(ikh, iks).para_bs.fix_mu = true;

                    self.xialpha(ikh, iks).para_alpha.no_zero_peak = false;
                    self.xialpha(ikh, iks).para_alpha.mu_df_lb = 0.5;
                    self.xialpha(ikh, iks).para_alpha.mu_df_ub = 3;

                    self.xialpha(ikh, iks).lambda = [0; 0; 1e-3; 1e-3; 2 * 1e-2]; %[1e-1;1e-1;1e-1]; %mu sigma startend no_negative
                    self.xialpha(ikh, iks).name = self.hos.name;
                    if isempty(self.verbose) && ispc %ispc false
                        self.verbose = true;
                    end

                    if self.verbose
                        self.xialpha(ikh, iks).verbose = self.verbose;
                    end

                end

            end

            self.xialpha = self.xialpha(:);
        end

        function fit(self)
            model = flip(self.xialpha);
            if ~isinparfor && length(model)>1
                startmatlabpool(min(feature('numcores'), self.maxworker));
                nidloop = length(model);
                idx = gen_loopgroup(length(model), nidloop);
                num_loop = numel(idx);

                for iloop = 1:num_loop
                    if self.verbose || nidloop == 1 %false %ispc %false %false ispc

                        for ik = idx{iloop}
                            disp(['[BiXiAlphaBootstrap]: ', num2str(ik), 'th fitting ', model(ik).hos_components, ' ks=', num2str(model(ik).ks), ' kh=', num2str(model(ik).kh)])
                            tempmodel = model(ik);
                            tempmodel.fit;
                            model(ik) = tempmodel;
                        end

                    else

                        parfor ik = idx{iloop}
                            disp(['[BiXiAlphaBootstrap parallel]: ', num2str(ik), 'th fitting ', model(ik).hos_components, ' ks=', num2str(model(ik).ks), ' kh=', num2str(model(ik).kh)])
                            tempmodel = model(ik);
                            tempmodel.fit;
                            model(ik) = tempmodel;
                        end

                    end

                end

            else

                for ik = 1:length(self.xialpha)
                    model(ik).fit;
                end

            end

            self.models = flip(model);
        end

        function pick_best(self)
            self.crt(:, 1) = cat_struct_fields(self.models, self.metric, false, true);

            if startsWith(self.metric, 'aic') || startsWith(self.metric, 'bic') || startsWith(self.metric, 'nllh') || startsWith(self.metric, 'mse')
                [~, idx_best] = min(self.crt, [], 1);
            else
                [~, idx_best] = max(self.crt, [], 1);
            end

            self.xialpha = copy(self.models(idx_best(:, 1)));
            disp(['best model is: ', ' ks=', num2str(self.xialpha.ks), ' kh=', num2str(self.xialpha.kh)])

        end

        function ttest_component(self, comp)
            disp(['chi-square test for ', comp])

            switch comp
                case 'raw'
                    bc = mean(self.xialpha.bc, 3);
                case 'fit'
                    bc = self.xialpha.bchat;
                case 'res'
                    bc = mean(self.xialpha.bc, 3) - self.xialpha.bchat;
                case 'xi'
                    bc = self.xialpha.bchat_component('xi');
                case 'alpha'
                    bc = self.xialpha.bchat_component('alpha');
                case 'xi_xi'
                    bc = self.xialpha.bchat_component('xi_xi');
                case 'alpha_alpha'
                    bc = self.xialpha.bchat_component('alpha_alpha');
                otherwise
                    error('wrong component')
            end

            top_level = 1;

            type_stats = {'median', 'max'};

            for itype = 1:length(type_stats)
                self.xialpha.stat.(comp).(type_stats{itype}) = calc_bc_stat(bc, type_stats{itype}, self.Ns_eff, self.p_threshold);

            end

        end

        function gen_null_nonpara(self, comp)
            % DEPRECATED: Bootstrap null generation — part of runfitnull pathway.
            % Retained for reference; not used in the main runfit() pathway.
            % The bootstrap null approach was abandoned because fitting 2*Nb
            % null models per channel is computationally prohibitive for
            % large-scale datasets (1772 iEEG channels, 960 EEG subjects).

            switch comp
                case 'xi'
                    self.null_xi.xialpha_null = copy(self.xialpha);
                    self.null_xi.xialpha_null.para_bs.kernel.para.h(1, :, :) = 0;
                    self.null_xi.xialpha_null.para_bs.kernel.para.h(:, 1, :) = 0;
                case 'alpha'
                    self.null_alpha.xialpha_null = copy(self.xialpha);
                    self.null_alpha.xialpha_null.para_bs.kernel.para.h(2:end, :, :) = 0;
                    self.null_alpha.xialpha_null.para_bs.kernel.para.h(:, 2:end, :) = 0;
                otherwise
                    error('wrong component')
            end

            field_comp = ['null_', comp];
            shat_null = self.(field_comp).xialpha_null.shat;
            bshat_null = self.(field_comp).xialpha_null.bshat;

            for ib = self.Nb:-1:1
                self.(field_comp).xialpha_sg(ib, :) = copy(self.xialpha);
                self.(field_comp).xialpha_sg(ib, :).name = append(self.(field_comp).xialpha_sg(ib, :).name, ['_', field_comp, '_', num2str(ib)]);
                self.(field_comp).xialpha_sg(ib, :).para_xi.kernel.para = self.xialpha.para_xi.kernel.para;
                self.(field_comp).xialpha_sg(ib, :).para_alpha.kernel.para = self.xialpha.para_alpha.kernel.para;
                self.(field_comp).xialpha_sg(ib, :).para_bs.kernel.para = self.xialpha.para_bs.kernel.para;

                % fix the kernel parameters except height h
                self.(field_comp).xialpha_sg(ib, :).para_xi.kernel.fix = self.xialpha.para_xi.kernel.para;
                self.(field_comp).xialpha_sg(ib, :).para_alpha.kernel.fix = self.xialpha.para_alpha.kernel.para;
                self.(field_comp).xialpha_sg(ib, :).para_bs.kernel.fix = self.xialpha.para_bs.kernel.para;

                self.(field_comp).xialpha_sg(ib, :).para_xi.kernel.fix.h = NaN(size(self.(field_comp).xialpha_sg(ib, :).para_xi.kernel.fix.h));
                self.(field_comp).xialpha_sg(ib, :).para_alpha.kernel.fix.h = NaN(size(self.(field_comp).xialpha_sg(ib, :).para_alpha.kernel.fix.h));
                self.(field_comp).xialpha_sg(ib, :).para_bs.kernel.fix.h = NaN(size(self.(field_comp).xialpha_sg(ib, :).para_bs.kernel.fix.h));
            end

            if self.hos.log_var_s
                sres = log(self.hos.s) - log(self.xialpha.shat);
            else
                sres = self.hos.s - self.xialpha.shat;
            end

            bsres = self.hos.bs - self.xialpha.bshat;

            for ib = 1:self.Nb
                self.(field_comp).idx_bootstrap(:, ib) = randi(self.hos.Ns, [self.Ns, 1]);
                idx_bootstrap = self.(field_comp).idx_bootstrap(:, ib); %#ok<*PROPLC>

                if self.multiple_seg

                    if self.hos.log_var_s
                        sres_star = sres(:, idx_bootstrap);
                        bsres_star = bsres(:, :, idx_bootstrap);
                        sres_star = sres_star - mean(sres_star, 2);
                        bsres_star = bsres_star - mean(bsres_star, 3);
                        sstar = exp(log(shat_null) + sres_star);
                        bsstar = bshat_null + tril2full(bsres_star, 0, false);
                        sstarvar = var(sres_star, [], 2);
                        bsstarvar = real(sum(bsres_star .* conj(bsres_star), 3) ./ self.hos.Ns);
                        bsdenom = self.hos.bsdenom(:, :, idx_bootstrap);
                    else
                        error('not implemented yet');
                    end

                else % fit single mean
                    warning('not fits to the bootstrap test for nonlinear regression, can not demean')

                    if self.hos.log_var_s
                        sstar = exp(log(shat_null) + mean(sres(:, idx_bootstrap)));
                        bsstar = bshat_null + tril2full(mean(bsres(:, :, idx_bootstrap), 0, false));
                        sstarvar = var(sres(:, idx_bootstrap), [], 2) ./ Ns; %#ok<*CPROPLC>
                        bsstarvar = real(sum(bsres(:, :, idx_bootstrap) .* conj(bsres(:, :, idx_bootstrap)), 3) ./ self.hos.Ns) ./ self.hos.Ns;
                    else
                        error('not implemented yet');
                    end

                end

                self.(field_comp).xialpha_sg(ib, :).hos.bsdenom = bsdenom; %set hos first, if not it will cover the svar and bsvar

                self.(field_comp).xialpha_sg(ib, :).s = sstar;
                self.(field_comp).xialpha_sg(ib, :).svar = sstarvar;

                self.(field_comp).xialpha_sg(ib, :).bs = bsstar;
                self.(field_comp).xialpha_sg(ib, :).bsvar = bsstarvar;
            end

        end

        function generate_null(self)
            % DEPRECATED: See gen_null_nonpara — part of runfitnull pathway.

            if self.verbose
                self.Nb = 1;
            end

            gen_null_nonpara(self, 'xi');
            gen_null_nonpara(self, 'alpha')
        end

        function fit_null(self)
            % DEPRECATED: Fits 2*Nb null models — part of runfitnull pathway.
            % Computationally prohibitive for large datasets.
            modelnull = [self.null_xi.xialpha_sg; self.null_alpha.xialpha_sg];

            if isinparfor || length(modelnull) == 1

                for ib = 1:(2 * self.Nb)
                    disp(['[BiXiAlphaBootstrap]: ', num2str(ib), 'th fitting ', modelnull(ib).hos_components])
                    modelnull(ib).fit;
                end

            else
                startmatlabpool(min(feature('numcores'), self.maxworker));
                nidloop = length(modelnull);
                idx = gen_loopgroup((2 * self.Nb), nidloop);
                num_loop = numel(idx);

                for iloop = 1:num_loop
                    parfor ib = idx{iloop}
                        disp(['[BiXiAlphaBootstrap parallel]: ', num2str(ib), 'th fitting ', modelnull(ib).hos_components])
                        tempmodel = modelnull(ib);
                        tempmodel.fit;
                        modelnull(ib) = tempmodel;
                    end

                end

            end

            self.modelsnull = modelnull;
        end

        function statistic_compare(self)
            % DEPRECATED: Bootstrap comparison — part of runfitnull pathway.
            self.null_xi.xialpha_sg = copy(self.modelsnull(1:self.Nb));
            self.null_alpha.xialpha_sg = copy(self.modelsnull(self.Nb + 1:end));

            bootstrap_test_component(self, 'xi');
            bootstrap_test_component(self, 'alpha');
            bootstrap_test_component(self, 'xi_xi');
            bootstrap_test_component(self, 'alpha_alpha');
        end

        function bootstrap_test_component(self, comp_test)
            % DEPRECATED: Per-component bootstrap test — part of runfitnull pathway.
            disp(['bootstrap test for ', comp_test])

            switch comp_test
                case 'xi'
                    comp = 'xi';
                    bc = self.xialpha.bchat_component('xi');
                case 'alpha'
                    comp = 'alpha';
                    bc = self.xialpha.bchat_component('alpha');
                case 'xi_xi'
                    comp = 'xi';
                    bc = self.xialpha.bchat_component('xi_xi');
                case 'alpha_alpha'
                    comp = 'alpha';
                    bc = self.xialpha.bchat_component('alpha_alpha');
                otherwise
                    error('wrong component')
            end

            field_comp = ['null_', comp];

            for ib = self.Nb:-1:1
                bc_null(:, :, ib) = abs(self.(field_comp).xialpha_sg(ib).bchat);
            end

            top_level = 1;
            type_stats = {'sum', 'mean', 'robustmean', 'logmean', 'robustlogmean', 'max', 'iqr'};

            for istat = 1:length(type_stats)
                self.(field_comp).stat.(comp_test).(type_stats{istat}) = calc_bc_stat(bc, type_stats{istat}, self.xialpha.hos.Ns, self.p_threshold, top_level);
            end

        end

        function stat_bc_compare(self, comp, stat_type, ispolar)
            % DEPRECATED: Bicoherence bootstrap comparison — part of runfitnull pathway.
            field_comp = ['null_', comp];
            for ib = self.Nb:-1:1

                switch stat_type
                    case 't'

                        if ispolar
                            bc_t_null(:, :, ib) = self.(field_comp).xialpha_sg(ib).stat.bc_polar.(comp).t;
                            bc_t = self.xialpha.stat.bc_polar.(comp).t;
                        else
                            bc_t_null(:, :, ib) = self.(field_comp).xialpha_sg(ib).stat.bc.(comp).t;
                            bc_t = self.xialpha.stat.bc.(comp).t;
                        end

                    case 'bc'

                        if ispolar
                            bc_t_null(:, :, ib) = abs(self.(field_comp).xialpha_sg(ib).bchat_component(comp));
                            bc_t = abs(self.xialpha.bchat_component(comp));
                        else
                            error('no such method')
                        end

                    otherwise
                        error('no such method')
                end

            end

            if ispolar
                bc_t_null = real(bc_t_null);
                bc_t = real(bc_t);
            else
                bc_t_null = cat(4, real(bc_t_null), imag(bc_t_null));
                bc_t = cat(4, real(bc_t), imag(bc_t));
            end

            if ispolar

                distribution = 'Nakagami';
            else
                distribution = 'Kernel';

            end

            self.(field_comp).bc.(comp) = self.calc_stat_parfor(bc_t, bc_t_null, 3, self.alpha_level, ispolar, distribution);
        end

        function stat = calc_stat_parfor(~, samples, surrogates, dim, alpha_level, singleside, varargin)
            % DEPRECATED: Parfor bootstrap statistic — part of runfitnull pathway.
            tsg = surrogates;
            t = samples;
            if ~exist('dim', 'var') || isempty(dim), dim = 1; end
            idx_dim_sg = 1:ndims(surrogates);
            idx_dim_sp = 1:ndims(samples);
            sz_sg = size(surrogates);
            sz_sp = size(samples);
            N = sz_sg(dim);
            surrogates = reshape(permute(surrogates, [idx_dim_sg(1:dim - 1), idx_dim_sg(dim + 1:end), idx_dim_sg(dim)]), [], sz_sg(dim));

            if dim <= ndims(samples)
                samples = reshape(permute(samples, [idx_dim_sp(1:dim - 1), idx_dim_sp(dim + 1:end), idx_dim_sp(dim)]), [], sz_sp(dim));
            else
                samples = reshape(samples, [], 1);
            end

            Nentry = size(surrogates, 1);
            para = varargin;

            pd_temp = fitdist(rand(10, 1), para{:});
            pd = repmat(get_empty_obj(pd_temp), Nentry, 1);
            th_low = NaN(Nentry, 1);
            th_up = NaN(Nentry, 1);
            p = NaN(Nentry, 1);
            h = NaN(Nentry, 1);

            parfor ientry = 1:Nentry
                surrogate = surrogates(ientry, :)';
                sample = samples(ientry, :);

                if any(isnan(surrogate)) ...
                        || any(isinf(surrogate)) ...
                        || all(surrogate == 0) ...
                        || any(isnan(sample)) ...
                        || any(isinf(sample))
                    continue;
                end

                if singleside
                    samplemin = min(sample(sample > 0));

                    if isempty(samplemin)
                        samplemin = inf;
                    end

                    surrogate(surrogate <= 0) = min(min(min(surrogate(surrogate > 0)), eps), samplemin);
                end

                ws = warning('off');
                pd(ientry, :) = fitdist(surrogate, para{:}); %#ok<PFBNS>
                warning(ws);
                th_low(ientry, :) = pd(ientry, :).icdf(alpha_level);
                th_up(ientry, :) = pd(ientry, :).icdf(1 - alpha_level);
                p(ientry, :) = pd(ientry, :).cdf(sample, 'upper');

                if ~singleside
                    p(ientry, :) = 2 * min(p(ientry, :), 1 - p(ientry, :));
                end

                h(ientry, :) = double(p(ientry, :) < alpha_level);
            end

            sz_stat = sz_sp;

            if dim <= ndims(samples)
                sz_stat(dim) = [];
                sz_stat = [sz_stat, sz_sp(dim)];
            else
                sz_stat = [sz_stat, 1];
            end

            idx_dim_stat = 1:length(sz_stat);
            idx_reorder = [idx_dim_stat(1:dim - 1), idx_dim_stat(end), idx_dim_stat(dim:end - 1)];
            pd = permute(reshape(pd, sz_stat), idx_reorder);
            th_low = permute(reshape(th_low, sz_stat), idx_reorder);
            th_up = permute(reshape(th_up, sz_stat), idx_reorder);
            p = permute(reshape(p, sz_stat), idx_reorder);
            h = permute(reshape(h, sz_stat), idx_reorder);
            h_all = any(h(:));

            stat = struct('pd', pd, 'tsg', tsg, 't', t, 'th_low', th_low, 'th_up', th_up, 'p', p, 'h', h, 'h_all', h_all);

        end

        function pd = fitdistmulti(~, x, dim, varargin)
            sz = size(x);
            if ~exist('dim', 'var') || isempty(dim), dim = 1; end
            sz(dim) = [];
            x = linspacemulti(ones(1, length(sz)), sz, sz);
            idx = get_ndgrid(x, 'list');
            ws = warning('off');

            for i = length(idx):-1:1
                pd(idx(i, :)) = fitdist(squeeze(x(idx(i, :))), varargin{:});
            end

            warning(ws);

        end

        function stat_para_compare(self, comp, ispolar)
            % DEPRECATED: Parameter bootstrap comparison — part of runfitnull pathway.
            field_comp = ['null_', comp];

            switch comp
                case {'xi'}
                    mask = false(size(self.xialpha.para_bs.kernel.para.h));
                    mask(1, 1, :) = true;
                case {'alpha'}
                    mask = true(size(self.xialpha.para_bs.kernel.para.h));
                    mask(1, :, :) = false;
                    mask(:, 1, :) = false;
            end

            self.(field_comp).para.tvec = self.xialpha.para_bs.kernel.para.h(mask);

            if ispolar
                Nh = round(length(self.(field_comp).para.tvec) / 2);
                self.(field_comp).para.tvec = sqrt(self.(field_comp).para.tvec(1:Nh) .^ 2 + self.(field_comp).para.tvec(Nh + 1:2 * Nh) .^ 2);
            end

            self.(field_comp).para.(comp).tvec = [];

            for ib = self.Nb:-1:1
                self.(field_comp).para.(comp).tvec(:, ib) = self.(field_comp).xialpha_sg(ib).para_bs.kernel.para.h(mask);

                if ispolar
                    tvec(:, ib) = sqrt(self.(field_comp).para.(comp).tvec(1:Nh, ib) .^ 2 + self.(field_comp).para.(comp).tvec(Nh + 1:2 * Nh, ib) .^ 2);
                end

            end

            self.(field_comp).para.(comp).tvec = tvec;

            Np = size(~isnan(sum(self.(field_comp).para.(comp).tvec, 2)), 1);

            self.(field_comp).para.(comp) = mv_rmfield(self.(field_comp).para.(comp), {'pd', 'th', 'p', 'h'});

            if ispolar
                distribution = 'Nakagami';
            else
                distribution = 'Kernel';
            end

            ws = warning('off');

            for ip = Np:-1:1
                tvec_sum = sum(self.(field_comp).para.(comp).tvec(ip, :));

                if isnan(tvec_sum) || tvec_sum == 0
                    continue
                end

                if ispolar
                    sample = self.(field_comp).para.tvec(ip, :);
                    samplemin = min(sample(sample > 0));

                    if isempty(samplemin)
                        samplemin = inf;
                    end

                    self.(field_comp).para.(comp).tvec(ip, self.(field_comp).para.(comp).tvec(ip, :) <= 0) = min(min(min(self.(field_comp).para.(comp).tvec(ip, self.(field_comp).para.(comp).tvec(ip, :) > 0)), eps), samplemin);
                end

                self.(field_comp).para.(comp).pd(ip, :) = fitdist(self.(field_comp).para.(comp).tvec(ip, :)', distribution); %,'Bandwidth', 0.5
                self.(field_comp).para.(comp).th(ip, 1) = self.(field_comp).para.(comp).pd(ip).icdf(self.alpha_level); %Support
                self.(field_comp).para.(comp).th(ip, 2) = self.(field_comp).para.(comp).pd(ip).icdf(1 - self.alpha_level);

                if ~ispolar
                    self.(field_comp).para.(comp).p(ip, :) = 2 * min(self.(field_comp).para.(comp).pd(ip).cdf(self.(field_comp).para.tvec(ip, :)), 1 - self.(field_comp).para.(comp).pd(ip).cdf(self.(field_comp).para.tvec(ip, :)));
                else
                    self.(field_comp).para.(comp).p(ip, :) = self.(field_comp).para.(comp).pd(ip).cdf(self.(field_comp).para.tvec(ip, :), 'upper');
                end

                self.(field_comp).para.(comp).hyp(ip, :) = self.(field_comp).para.(comp).p(ip, :) < self.alpha_level;
            end

            warning(ws);

            if length(self.(field_comp).para.(comp).pd) < Np
                self.(field_comp).para.(comp).pd(Np, :) = get_empty_obj(self.(field_comp).para.(comp).pd(1));
                self.(field_comp).para.(comp).th(Np, 1) = get_empty_obj(self.(field_comp).para.(comp).th(1, 1));
                self.(field_comp).para.(comp).th(Np, 2) = get_empty_obj(self.(field_comp).para.(comp).th(1, 2));
                self.(field_comp).para.(comp).p(Np, :) = get_empty_obj(self.(field_comp).para.(comp).p(1, :));
                self.(field_comp).para.(comp).hyp(Np, :) = get_empty_obj(self.(field_comp).para.(comp).hyp(1, :));
            end

            if ispolar
                mask = mask(:, :, 1);
            end

            self.(field_comp).para.(comp).pd_mat = bmap(self.(field_comp).para.(comp).pd, mask, mask);
            self.(field_comp).para.(comp).th_mat_low = bmap(self.(field_comp).para.(comp).th(:, 1), mask, mask);
            self.(field_comp).para.(comp).th_mat_up = bmap(self.(field_comp).para.(comp).th(:, 2), mask, mask);
            self.(field_comp).para.(comp).p_mat = bmap(self.(field_comp).para.(comp).p, mask, mask);
            self.(field_comp).para.(comp).hyp_mat = bmap(self.(field_comp).para.(comp).hyp, mask, mask);
            self.(field_comp).para.tmat = bmap(self.(field_comp).para.tvec, mask, mask);
        end

        function compact_result(self)
            self.hos = [];
            self.models = [];
            self.modelsnull = [];
            self.xialpha.temp = [];
            self.xialpha.model = [];

            if ~isempty(self.null_xi)
                self.null_xi.xialpha_null = [];

                self.null_xi.xialpha_sg(1:end) = []; %keep one null

            end

            if ~isempty(self.null_alpha)
                self.null_alpha.xialpha_null = [];

                self.null_alpha.xialpha_sg(1:end) = []; %keep one null
            end

        end

        function runfit(self)

            if isempty(self.rngseed)
                self.rngseed = randi(1000000);
            end

            rng(self.rngseed);

            self.init();
            tic;
            self.fit();

            self.pick_best();
            if isempty(self.Ns_eff)
                if isempty(self.x)
                    error('needs the data x for estimate Ns_eff');
                end
                self.Ns_eff = get_hos_Ns_eff(self.x, self.xialpha.hos, 100, false, true);
            end
            self.ttest_component('raw');
            self.ttest_component('xi');
            self.ttest_component('alpha');


            self.logging.t = toc;
            disp(['[BiXiAlphaBootstrap] Fit orignal data elapsed time=', num2str(self.logging.t)])
        end

        function runfitnull(self)
            % DEPRECATED: Bootstrap null testing pathway.
            % Abandoned because fitting 2*Nb null models per channel is
            % computationally prohibitive for large datasets (1772 iEEG
            % channels, 960 EEG subjects).
            % The main statistical testing pathway is runfit() ->
            % ttest_component(), which uses asymptotic chi-squared tests.
            self.generate_null();

            if self.tf_runtest %~ispc ~false
                tic;
                self.fit_null();
                self.logging.t = toc;
                disp(['[BiXiAlphaBootstrap] Fit null data elapsed time=', num2str(self.logging.t)])
                self.statistic_compare();
            end

        end

        function runall(self)

            self.runfit();
            self.runfitnull();
            self.compact_result();
        end

        function show(self)
            self.xialpha.show;
        end

        function show_stat(self)
            comp = {'xi', 'alpha'};

            for icomp = 1:length(comp)
                subplot(2, 1, icomp)
                hold on;
                field_comp = ['null_', comp{icomp}];

                if ~isfield(self.(field_comp), (comp{icomp}))
                    continue;
                end

                Np = size(self.(field_comp).(comp{icomp}).pd, 1);
                x = linspace(min(self.(field_comp).para.tvec, [], 'all'), max(self.(field_comp).para.tvec, [], 'all'), 100);
                cm = colormap(lines(Np));

                for ip = 1:Np

                    if isnan(sum(self.(field_comp).para.tvec(ip, :))) || isempty(self.(field_comp).(comp{icomp}).pd(ip))
                        continue
                    end

                    self.(field_comp).(comp{icomp}).pd(ip, :).plot;
                    xline(self.(field_comp).para.tvec(ip, :), 'Color', cm(ip, :), 'LineWidth', 2)
                end

            end

        end

        function show_stat_matrix(self, comp)
            field_comp = ['null_', comp];

            if ~isfield(self.(field_comp).para, (comp))
                return;
            end

            clf
            mu1 = [0, get_harmonic_mu_sigma(self.xialpha.para_bs.kernel.para.mu(end), 0, self.xialpha.ks, self.xialpha.kh, self.xialpha.para_alpha.no_zero_peak)];
            [mux, muy] = meshgrid(mu1, mu1);

            ispolar = size(self.(field_comp).para.(comp).pd_mat, 3) == 1;

            if strcmpi(comp, 'xi')

                if ~ispolar
                    subplot(2, 1, 1)
                    plot_para_dist(self, self.(field_comp).para.(comp).pd_mat(1, 1, 1), self.(field_comp).para.tmat(1, 1, 1), self.(field_comp).para.(comp).hyp_mat(1, 1, 1)); %self.xialpha.(field_comp).para.tvec(1)
                    subplot(2, 1, 2)
                    plot_para_dist(self, self.(field_comp).para.(comp).pd_mat(1, 1, 2), self.(field_comp).para.tmat(1, 1, 2), self.(field_comp).para.(comp).hyp_mat(1, 1, 2));
                else
                    plot_para_dist(self, self.(field_comp).para.(comp).pd_mat(1, 1, 1), self.(field_comp).para.tmat(1, 1, 1), self.(field_comp).para.(comp).hyp_mat(1, 1, 1));
                end

            elseif strcmpi(comp, 'alpha')

                for k = 1:size(self.(field_comp).para.(comp).pd_mat, 3)

                    for j = 1:size(self.(field_comp).para.(comp).pd_mat, 2)

                        for i = 1:size(self.(field_comp).para.(comp).pd_mat, 1)

                            if ~ispolar
                                isubplot = 2 * (i - 1) * size(self.(field_comp).para.(comp).pd_mat, 2) + j + (k - 1) * size(self.(field_comp).para.(comp).pd_mat, 2);
                                subplot(size(self.(field_comp).para.(comp).pd_mat, 1), size(self.(field_comp).para.(comp).pd_mat, 2) * 2, isubplot)
                            else
                                isubplot = (i - 1) * size(self.(field_comp).para.(comp).pd_mat, 2) + j + (k - 1) * size(self.(field_comp).para.(comp).pd_mat, 2);
                                subplot(size(self.(field_comp).para.(comp).pd_mat, 1), size(self.(field_comp).para.(comp).pd_mat, 2), isubplot)
                            end

                            plot_para_dist(self, self.(field_comp).para.(comp).pd_mat(i, j, k), self.(field_comp).para.tmat(i, j, k), self.(field_comp).para.(comp).hyp_mat(i, j, k)); %self.xialpha.(field_comp).para.tmat(i, j, k)
                            title([num2str(mux(i, j), '%.1f'), ',', num2str(muy(i, j), '%.1f')])
                        end

                    end

                end

            end

        end

        function plot_para_dist(~, pd, tval, hval)

            if strcmpi(pd.DistributionName, 'Uniform')
                grid off
                box off
                return
            end

            hold on;

            if hval
                h = pd.plot;
                h(2).FaceColor = 'red';
                xline(tval, 'LineWidth', 2, 'Color', 'blue')

            else
                h = pd.plot;
                h(2).FaceColor = 'green';
                xline(tval, 'LineWidth', 2, 'Color', 'blue')
            end

            hold off;
        end

        function disp_stat(self)
            comps = fieldnames(self.xialpha.stat);

            for icomp = 1:length(comps)
                comp = comps{icomp};
                disp(['chi-square test for ', comp])
                type_stats = fieldnames(self.xialpha.stat.(comp));

                for istat = 1:length(type_stats)
                    stat = self.xialpha.stat.(comp).(type_stats{istat});
                    type_stat = type_stats{istat};

                    if stat.h
                        disp(['Test for ', type_stat, ' is significant', ' p=', num2str(stat.p), ', this process is "needed"'])
                    else
                        disp(['Test for ', type_stat, ' is not significant', ' p=', num2str(stat.p), ', this process is "not needed"'])
                    end

                end

            end

            comps_test = {'xi', 'alpha', 'xi_xi', 'alpha_alpha'};

            for icomp = 1:length(comps_test)
                comp_test = comps_test{icomp};

                switch comp_test
                    case {'xi', 'xi_xi'}
                        comp = 'xi';
                    case {'alpha', 'alpha_alpha'}
                        comp = 'alpha';
                end

                field_comp = ['null_', comp];

                disp(['bootstrap test for ', comp_test])

                for istat = 1:length(type_stats)

                    if isempty(self.(field_comp))
                        disp('run bootstrap first');
                        continue;
                    end

                    stat = self.(field_comp).stat.(comp_test).(type_stats{istat});
                    type_stat = type_stats{istat};

                    if stat.h
                        disp(['Test for ', type_stat, ' is significant', ' p=', num2str(stat.p), ', this process is "needed"'])
                    else
                        disp(['Test for ', type_stat, ' is not significant', ' p=', num2str(stat.p), ', this process is "not needed"'])
                    end

                end

            end

        end

        function plot_null_hypo_dist(self)
            hFig = gcf;
            hFig.WindowState = 'maximized';

            if isfield(self.null_xi, 'xi') && isfield(self.null_xi.para.xi, 'pd') && ~isempty(self.null_xi.para.xi.pd)
                subplot(3, 2, 1)
                pdfplot(self.null_xi.para.xi.pd);
                xline(self.null_xi.para.bchat_xi_max)
                title('xi max distribution')
            end

            if isfield(self.null_alpha.para.alpha, 'pd') && ~isempty(self.null_alpha.para.alpha.pd)
                subplot(3, 2, 2)
                pdfplot(self.null_alpha.para.alpha.pd);
                xline(self.null_alpha.para.bchat_alpha_max)
                title('alpha max distribution')
            end

        end

        function summary = get.summary(self)

            if length(self.xialpha) == 1
                summary = [self.xialpha.summary];
            else
                summary = [];
            end

            if isfield(self.null_xi, 'xi_xi') && isfield(self.null_xi.xi_xi, 'h') && ~isempty(self.null_xi.xi_xi.h)
                stat.xi_h = self.null_xi.xi.h;
                stat.xi_xi_h = self.null_xi.xi_xi.h;
                stat.xi_alpha_h = self.null_xi.xi_alpha.h;

                stat.xi_p = self.null_xi.xi.p;
                stat.xi_xi_p = self.null_xi.xi_xi.p;
                stat.xi_alpha_p = self.null_xi.xi_alpha.p;

                stat.bchat_xi_max = self.null_xi.bchat_xi_max;
                stat.bchat_xi_xi_max = self.null_xi.bchat_xi_xi_max;
                stat.bchat_xi_alpha_max = self.null_xi.bchat_xi_alpha_max;

                stat.null_bchat_xi_max = max(self.null_xi.xi.bchatmax);
                stat.null_bchat_xi_xi_max = max(self.null_xi.xi_xi.bchatmax);
                stat.null_bchat_xi_alpha_max = max(self.null_xi.xi_alpha.bchatmax);

                summary = [summary, struct2table(stat)];
            end

        end

    end

end