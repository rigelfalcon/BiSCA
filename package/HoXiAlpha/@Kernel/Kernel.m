classdef Kernel %< handle & matlab.mixin.Copyable

    properties
        kfun
        type
        para % a struct of parameters
        lb
        ub
        fix
        t
        info
        use_info = false;
        hist
    end

    properties (Hidden)
        vpara % a vector of parameters
        vlb % a vector of lower bounds of parameters
        vub
        vfix
        vt
        para_info
        lb_info
        ub_info
        fix_info
        t_info
    end

    methods

        function self = Kernel(type, para)

            if nargin == 0
                self.type = 'gaussian_kernel';
                self.para = struct('h', 1, 'mu', 0, 'sigma', 1);
            elseif nargin == 1
                self.type = type;
                self = self.init_para(type);
            elseif nargin == 2
                self.type = type;
                self.para = para;
            end

        end

        %% set methods
        function self = set.para(self, para)
            self.para = para;
            [~, self.para_info] = struct2vec(para);
        end

        function self = set.lb(self, lb)
            self.lb = lb;
            [~, self.lb_info] = struct2vec(lb);
        end

        function self = set.ub(self, ub)
            self.ub = ub;
            [~, self.ub_info] = struct2vec(ub);
        end

        function self = set.t(self, t)
            self.t = t;
            [~, self.t_info] = struct2vec(t);
        end

        function self = set.vpara(self, value)
            self.vpara = value;
            % set all parameters from a vector
            self.para = vec2struct(value, self.para_info);
        end

        function self = set.vlb(self, value)
            self.vlb = value;

            if isempty(self.lb_info)
                self.lb_info = self.para_info;
            end

            % set all lower bounds of parameters from a vector
            self.lb = vec2struct(value, self.lb_info);
        end

        function self = set.vub(self, value)
            self.vub = value;

            if isempty(self.ub_info)
                self.ub_info = self.para_info;
            end

            % set all upper bounds of parameters from a vector
            self.ub = vec2struct(value, self.ub_info);
        end

        function self = set.vfix(self, value)
            self.vfix = value;

            if isempty(self.fix_info)
                self.fix_info = self.para_info;
            end

            % set all fix value of parameters from a vector
            self.fix = vec2struct(value, self.fix_info);
        end

        function self = set.vt(self, value)
            self.vt = value;

            if isempty(self.t_info)
                self.t_info = self.para_info;
            end

            % set all t value of parameters from a vector
            self.t = vec2struct(value, self.t_info);
        end

        function self = set.type(self, value)
            self.type = value;
            self.kfun = str2func(['Kernel.', value]); %#ok<*MCSUP>
        end

        %% get methods
        function vpara = get.vpara(self)
            % concatenate all parameters to a vector
            [vpara] = struct2vec(self.para); %,self.para_info
        end

        function vlb = get.vlb(self)
            % concatenate all lower bounds of parameters to a vector
            [vlb] = struct2vec(self.lb); %,self.lb_info
        end

        function vub = get.vub(self)
            % concatenate all upper bounds of parameters to a vector
            [vub] = struct2vec(self.ub); %,self.ub_info
        end

        function vfix = get.vfix(self)

            if ~isempty(self.fix)
                % concatenate all fix value of parameters to a vector
                [vfix] = struct2vec(self.fix); %,self.fix_info
            end

        end

        function vt = get.vt(self)
            % concatenate all t value of parameters to a vector
            [vt] = struct2vec(self.t); %,self.t_info
        end

        %% other methods
        function [y, ycomp] = eval(self, x)
            % feed structure para as cell format to kfun to get y
            para = struct2cell(self.para); %#ok<*PROPLC>
            if self.use_info || strcmpi(self.type, 'harmonic_t_flattop')
                info = struct2cell(self.info);
            end
            if nargout == 1
                if self.use_info || strcmpi(self.type, 'harmonic_t_flattop')
                    y = self.kfun(x, para{:}, info{:});
                else
                    y = self.kfun(x, para{:});
                end
            else
                if self.use_info
                    [y, ycomp] = self.kfun(x, para{:}, info{:});
                else
                    [y, ycomp] = self.kfun(x, para{:});
                end
            end

        end

        function T = check(self)
            T = [self.lb, self.para, self.ub];
            % check if the parameters are within the bounds
            % 0: not within bounds, 1: within bounds

            if ~isempty(self.vlb) && ~isempty(self.vpara) && ~isempty(self.vub)
                paras = [self.vlb, self.vpara, self.vub];
                tf = [paras(:, 1) <= paras(:, 2), paras(:, 2) <= paras(:, 3)];

                if any(tf(:) == 0)
                    disp('out of bounds')
                else
                    disp('within bounds')
                end

            end

        end

        function self = init_para(self, type)

            switch type
                case 'gaussian_kernel'
                    self.para = struct('h', 1, 'mu', 0, 'sigma', 1);
                case 'lorentzian_kernel'
                    self.para = struct('h', 1, 'mu', 0, 'k', 1, 'chi', 1);
                case 'linear_kernel'
                    self.para = struct('h', 1);
                case 'polynomial_kernel'
                    self.para = struct('h', 1, 'd', 2);
                case 'exponential_kernel'
                    self.para = struct('h', 1, 'beta', 1);
                case 'studentt_kernel'
                    self.para = struct('A', 1, 'B', 1, 'nu', 1, 'U', 0);
                case 'laplacian_kernel'
                    self.para = struct('h', 1, 'beta', 1);
                case 'periodic_kernel'
                    self.para = struct('h', 1, 'beta', 1, 'period', 1);
                case 'matern_kernel'
                    self.para = struct('h', 1, 'beta', 1, 'nu', 1);
                case 'harmonic_t_flattop'
                    self.use_info = true;
                    self.para = struct('h', [1, 0; 0, 2], 'mu', [0; 10], 'sigma', [10, 1; 1, 1], 'nu', [1; 1], 'd', [1, 1], 'b', [1, 1]);
                    self.info = struct('ks', 0, 'kh', 1, 'idx_sigma', [1, 0; 0, 2]);
                case 'harmonic_t_free'
                    self.use_info = true;
                    self.para = struct('h', [1, 0; 0, 2], 'mu', [0; 10], 'sigma', [10, 1; 1, 1], 'nu', [1; 1], 'd', [1, 1], 'b', [1, 1]);
                    self.info = struct('ks', 0, 'kh', 1, 'idx_sigma', [1, 0; 0, 2]);
                otherwise
                    error('no such method')
            end

        end

    end

    methods (Static)

        function para = setVectorField(para, value)
            % set all parameters of fields from a vector
            % 1. get the field names
            fnames = fieldnames(para);
            % 2. get the length of each field
            flength = cellfun(@(x) length(para.(x)), fnames);
            % 3. set the value of each field from the vector with vectorization
            for i = 1:length(flength)

                if i == 1
                    para.(fnames{i}) = value(1:flength(i));
                else
                    para.(fnames{i}) = value(sum(flength(1:i - 1)) + 1:sum(flength(1:i)));
                end

            end

        end

        function [y, ycomp] = gaussian_kernel(x, h, mu, sigma)
            if nargout == 1
                y = sum(h' .* exp(- (x - mu') .^ 2 ./ (2 * sigma' .^ 2)), 2);
            elseif nargout == 2
                ycomp = h' .* exp(- (x - mu') .^ 2 ./ (2 * sigma' .^ 2));
                y = sum(ycomp, 2);
            end
        end

        function [y, ycomp] = lorentzian_kernel(x, h, mu, k, chi)
            if nargout == 1
                y = sum(h' ./ (k + (x - mu') .^ (chi')), 2);
            elseif nargout == 2
                ycomp = h' ./ (k + (x - mu') .^ (chi'));
                y = sum(ycomp, 2);
            end
        end

        function y = linear_kernel(x, h)
            y = h * x;
        end

        function y = polynomial_kernel(x, h, d)
            y = h * x .^ d;
        end

        function y = exponential_kernel(x, h, beta)
            y = h * exp(-beta * x);
        end

        function y = studentt_kernel(x, A, B, nu, U)
            % follow the definition of student-t kernel in the Pascual's paper
            y = A ./ ((1.0 + ((x - U) ./ B) .^ 2) .^ nu);
        end

        function y = laplacian_kernel(x, h, beta)
            y = h * exp(-beta * abs(x));
        end

        function y = periodic_kernel(x, h, beta, period)
            y = h * exp(-2 * beta * sin(pi * abs(x) / period) .^ 2);
        end

        function y = matern_kernel(x, h, beta, nu)
            y = h * (1 + sqrt(2 * nu) / beta * abs(x)) .* exp(-sqrt(2 * nu) / beta * abs(x));
        end

        function [y, ycomp] = harmonic_t_flattop(x, h, mu, sigma, nu, d, b, ks, kh, idx_sigma, no_zero_peak, df, taperinfo)

            mu_in = mu;
            sigma_in = sigma;
            h_in = h;
            idx_sigma_in = idx_sigma;

            if nargin < 8 || isempty(idx_sigma)
                idx_sigma = [];
            end

            if nargin < 9 || isempty(no_zero_peak)
                no_zero_peak = false;
            end

            Nk = size(mu, 1); % number of different centers of kernels

            if Nk < length(h)

                for ik = Nk:-1:1

                    if mu(ik) ~= 0
                        mus{ik} = get_harmonic_mu_sigma(mu(ik), 0, ks, kh, no_zero_peak);
                    else
                        mus{ik} = mu(ik);
                    end

                end

                mu = cell2mat(mus)';
            end

            Nd = size(sigma, 2); % dimension of the input data

            Nslice = size(h, Nd + 1); % number of slices of the kernel, used for real and imaginary part of the kernel

            if ndims(x) ~= Nd && size(x, 2) ~= 1
                Nx = size(x);
                Nx = Nx(1:end - 1);
                x = reshape(x, [], Nd);
            else
                Nx = [];
            end

            N = size(x, 1);

            % remove zero height kernel
            if Nd == 2
                h = permute(h, [2, 1, 3]);
                [mux, muy] = meshgrid(mu, mu);
                mask = logical(sum(abs(h) ~= 0, 3));

                mux = mux(mask);
                muy = muy(mask);

                found = true(size(mux, 1), 1);
                mu = [mux, muy];
                mu(~found, :) = [];

                if ~isempty(idx_sigma)
                    idx_sigma = idx_sigma(mask);
                    idx_sigma(~found) = [];
                end

                h = h(logical(mask .* ones([ones(1, Nd), Nslice])));
                h = reshape(h, [], Nslice);
                h(~found, :) = [];
            end

            % remove kernel outside the region

            idx_outside = sum(mu, 2) > max(sum(x, 2));
            if size(sigma, 1) == size(mu, 1)
                sigma(idx_outside, :) = [];
            end

            if size(nu, 1) == size(mu, 1)
                nu(idx_outside, :) = [];
            end

            if size(d, 1) == size(mu, 1)
                d(idx_outside, :) = [];
            end

            if size(b, 1) == size(mu, 1)
                b(idx_outside, :) = [];
            end

            mu(idx_outside, :) = [];
            h(idx_outside, :) = [];

            if ~isempty(idx_sigma) && any(idx_outside)
                idx_sigma(idx_outside, :) = [];
            end

            if isempty(h)
                sz = size(x);
                ycomp = zeros([sz(1:(Nd - 1)), Nk]);
                y = zeros([sz(1:(Nd - 1)), 1]);
                return
            end

            h = permute(h, [3, 1, 2]);
            Sigma = sigma;
            Sigma(:, 2:end) = sigma(:, 1) .* sigma(:, 2:end);

            if Nd == 2
                Sigma = [Sigma(:, 1), Sigma];

                % Nu/d/b can be provided as either N-by-1 (shared across dims)
                % or N-by-2 (separate for f1/f2). Expand to N-by-3 for (f1,f2,f1+f2).
                if size(nu, 2) == 1
                    nu = repmat(nu, 1, 3);
                else
                    nu = [nu, mean(nu, 2)];
                end

                if size(d, 2) == 1
                    d = repmat(d, 1, 3);
                else
                    d = [d, mean(d, 2)];
                end

                if size(b, 2) == 1
                    b = repmat(b, 1, 3);
                else
                    b = [b, mean(b, 2)];
                end
            end

            if ~isempty(idx_sigma)
                Sigma = Sigma(idx_sigma, :);
                nu = nu(idx_sigma, :);
                d = d(idx_sigma, :);
                b = b(idx_sigma, :);
            end

            if nargout > 1
                [y, ycomp] = kernel_harmonic_t_flattop(x, h, mu, Sigma, nu, d, b);
            else
                [y] = kernel_harmonic_t_flattop(x, h, mu, Sigma, nu, d, b);
            end

            if size(y, 3) > 1
                y = y(:, :, 1) + 1i .* y(:, :, 2);
            end

            if ~isempty(Nx)
                y = reshape(y, [Nx(:)', 1]);
                y = tril2full(y, false, false);

            end

        end

        function [y, ycomp] = harmonic_t_free(x, h, mu, sigma, nu, d, b, ks, kh, idx_sigma, no_zero_peak, df, taperinfo)

            Nd = size(sigma, 2); % dimension of the input data

            Nslice = size(h, Nd + 1); % number of slices of the kernel, used for real and imaginary part of the kernel

            if ndims(x) ~= Nd && size(x, 2) ~= 1
                Nx = size(x);
                Nx = Nx(1:end - 1);
                x = reshape(x, [], Nd);
            else
                Nx = [];
            end

            N = size(x, 1);

            % remove zero height kernel
            if Nd == 2
                h = reshape(h, [], 2);
                idx_nzero = any(h ~= 0, 2);
                h = h(idx_nzero, :);
                mu = mu(idx_nzero, :);
                sigma = sigma(idx_nzero, :);
                nu = nu(idx_nzero, :);
                d = d(idx_nzero, :);
                b = b(idx_nzero, :);
            end

            if isempty(h)
                sz = size(x);
                ycomp = zeros([sz(1:(Nd - 1)), Nk]);
                y = zeros([sz(1:(Nd - 1)), 1]);
                return
            end

            if Nd == 2
                sigma = [sigma, mean(sigma, 2)];
                nu = [nu, mean(nu, 2)];
                d = [d, mean(d, 2)];
            end

            h = permute(h, [3, 1, 2]);

            if nargout > 1
                [y, ycomp] = kernel_harmonic_t_flattop(x, h, mu, sigma, nu, d, b);
            else
                [y] = kernel_harmonic_t_flattop(x, h, mu, sigma, nu, d, b);
            end

            if size(y, 3) > 1
                y = y(:, :, 1) + 1i .* y(:, :, 2);
            end

            if ~isempty(Nx)
                y = reshape(y, [Nx(:)', 1]);
                y = tril2full(y, false, false);

            end

        end
    end

end