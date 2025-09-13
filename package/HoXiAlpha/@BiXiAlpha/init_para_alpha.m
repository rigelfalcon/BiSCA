function init_para_alpha(self)

    switch self.para_alpha.type
        case 'tstudent'

            if strcmpi(self.para_alpha.type, 'tstudent') && strcmpi(self.peak_relation, 'harmonic_old')

                if ~isfield(self.para_alpha, 'kernel') || isempty(self.para_alpha.kernel) || ~strcmpi(self.para_alpha.kernel.type, self.kernel_type) || isempty(self.para_alpha.kernel.para) || length(self.para_alpha.kernel.para.h) ~= self.kmax
                    self.para_alpha.kernel = Kernel(self.kernel_type);
                    self.para_alpha.kernel.info.no_zero_peak = self.para_alpha.no_zero_peak;
                    self.k0 = double(~self.para_alpha.kernel.info.no_zero_peak);

                    self.para_alpha.kernel.para.h = [ones(1, self.kmax)];
                    self.para_alpha.kernel.para.mu = 11;
                    self.para_alpha.kernel.para.sigma = 1.39; %10 %
                    self.para_alpha.kernel.para.nu = 2; % 2 %60
                    self.para_alpha.kernel.para.d = 6.14; %3
                    self.para_alpha.kernel.para.b = 1;
                    self.para_alpha.kernel.info.ks = self.ks;
                    self.para_alpha.kernel.info.kh = self.kh;
                    self.para_alpha.kernel.info.idx_sigma = [];
                    self.para_alpha.kernel.info.df = self.df;
                    self.para_alpha.kernel.info.taperinfo = self.taperinfo;
                end

            elseif strcmpi(self.para_alpha.type, 'tstudent') && strcmpi(self.peak_relation, 'free')

                if ~isfield(self.para_alpha, 'kernel') || isempty(self.para_alpha.kernel) || ~strcmpi(self.para_alpha.kernel.type, self.kernel_type) || isempty(self.para_alpha.kernel.para) || length(self.para_alpha.kernel.para.h) ~= self.kmax
                    self.para_alpha.kernel = Kernel(self.kernel_type);
                    self.para_alpha.kernel.info.no_zero_peak = self.para_alpha.no_zero_peak;
                    self.k0 = double(~self.para_alpha.kernel.info.no_zero_peak);

                    self.para_alpha.kernel.para.h = [ones(1, self.kmax)];
                    self.para_alpha.kernel.para.mu = 11 .* ones(self.kmax, 1);
                    self.para_alpha.mu_init = self.para_alpha.kernel.para.mu;
                    if isfield(self.para_alpha, 'perpeak_shape') && self.para_alpha.perpeak_shape
                        shape_rep = ones(self.kmax, 1);
                    else
                        shape_rep = 1;
                    end
                    self.para_alpha.kernel.para.sigma = 1.39 .* shape_rep; %10 %
                    self.para_alpha.kernel.para.nu = 2 .* shape_rep; % 2 %60
                    self.para_alpha.kernel.para.d = 6.14 .* shape_rep; %3
                    self.para_alpha.kernel.para.b = 1 .* shape_rep;
                    self.para_alpha.kernel.info.ks = self.ks;
                    self.para_alpha.kernel.info.kh = self.kh;
                    self.para_alpha.kernel.info.idx_sigma = [];
                    self.para_alpha.kernel.info.df = self.df;
                    self.para_alpha.kernel.info.taperinfo = self.taperinfo;
                end
            elseif strcmpi(self.para_alpha.type, 'tstudent') && strcmpi(self.peak_relation, 'harmonic')

                if ~isfield(self.para_alpha, 'kernel') || isempty(self.para_alpha.kernel) || ~strcmpi(self.para_alpha.kernel.type, self.kernel_type) || isempty(self.para_alpha.kernel.para) || length(self.para_alpha.kernel.para.h) ~= self.kmax
                    self.para_alpha.kernel = Kernel(self.kernel_type);
                    self.para_alpha.kernel.info.no_zero_peak = self.para_alpha.no_zero_peak;
                    self.k0 = double(~self.para_alpha.kernel.info.no_zero_peak);

                    self.para_alpha.kernel.para.h = [ones(1, self.kmax)];
                    self.para_alpha.kernel.para.mu = 11 .* ones(self.kmax, 1);
                    self.para_alpha.mu_init = self.para_alpha.kernel.para.mu;
                    if isfield(self.para_alpha, 'perpeak_shape') && self.para_alpha.perpeak_shape
                        shape_rep = ones(self.kmax, 1);
                    else
                        shape_rep = 1;
                    end
                    self.para_alpha.kernel.para.sigma = 1.39 .* shape_rep; %10 %
                    self.para_alpha.kernel.para.nu = 2 .* shape_rep; % 2 %60
                    self.para_alpha.kernel.para.d = 6.14 .* shape_rep; %3
                    self.para_alpha.kernel.para.b = 1 .* shape_rep;
                    self.para_alpha.kernel.info.ks = self.ks;
                    self.para_alpha.kernel.info.kh = self.kh;
                    self.para_alpha.kernel.info.idx_sigma = [];
                    self.para_alpha.kernel.info.df = self.df;
                    self.para_alpha.kernel.info.taperinfo = self.taperinfo;
                end
            else
                error('no such peak_relation')
            end

        case 'gaussian'

            if ~isfield(self.para_alpha, 'kernel') || isempty(self.para_alpha.kernel) || ~strcmpi(self.para_alpha.kernel.type, self.kernel_type) || isempty(self.para_alpha.kernel.para) || length(self.para_alpha.kernel.para.h) ~= self.kmax
                self.para_alpha.kernel = Kernel('gaussian_kernel');
                self.para_alpha.kernel.info.no_zero_peak = self.para_alpha.no_zero_peak;
                self.k0 = double(~self.para_alpha.kernel.info.no_zero_peak);

                self.para_alpha.kernel.para.h = [ones(1, self.kmax)];
                self.para_alpha.kernel.para.mu = 11 .* ones(self.kmax, 1);
                self.para_alpha.mu_init = self.para_alpha.kernel.para.mu;
                self.para_alpha.kernel.para.sigma = 0.4 .* ones(self.kmax, 1); %10 %
                self.para_alpha.kernel.info.ks = self.ks;
                self.para_alpha.kernel.info.kh = self.kh;
                self.para_alpha.kernel.info.idx_sigma = [];
                self.para_alpha.kernel.info.df = self.df;
                self.para_alpha.kernel.info.taperinfo = self.taperinfo;
            end

        otherwise
            error('no such method')
    end

end