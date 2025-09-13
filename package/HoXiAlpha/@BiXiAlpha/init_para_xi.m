function init_para_xi(self)

self.para_xi.f_max = 3;
idx_xi_range = self.f < self.para_xi.f_max;

switch self.para_xi.type
    case 'tstudent'

        if ~isfield(self.para_xi, 'kernel') || isempty(self.para_xi.kernel) || ~strcmpi(self.para_xi.kernel.type, self.kernel_type) || isempty(self.para_xi.kernel.para) || length(self.para_xi.kernel.para.h) ~= 1
            self.para_xi.kernel = Kernel(self.kernel_type);
            self.para_xi.kernel.para.h = max(self.s, [], 'all');
            self.para_xi.kernel.para.mu = 0;
            self.para_xi.kernel.para.sigma = 1.70; %50; %5.^2
            self.para_xi.kernel.para.nu = 0.1; %3.2; %1
            self.para_xi.kernel.para.d = 27.18; %3;
            self.para_xi.kernel.para.b = 1; %3;
            self.para_xi.kernel.info.ks = 0;
            self.para_xi.kernel.info.kh = 1;
            self.para_xi.kernel.info.idx_sigma = [];
            self.para_xi.kernel.info.no_zero_peak = true;
            self.para_xi.kernel.info.df = self.df;
            self.para_xi.kernel.info.taperinfo = self.taperinfo;
        end
    case 'lorentzian'

        if ~isfield(self.para_xi, 'kernel') || isempty(self.para_xi.kernel) || ~strcmpi(self.para_xi.kernel.type, self.kernel_type) || isempty(self.para_xi.kernel.para) || length(self.para_xi.kernel.para.h) ~= 1
            self.para_xi.kernel = Kernel('lorentzian_kernel');
            self.para_xi.kernel.para.h = max(self.s, [], 'all');
            self.para_xi.kernel.para.mu = 0;
            self.para_xi.kernel.para.k = 1; %50; %5.^2
            self.para_xi.kernel.para.chi = 1; %3.2; %1
            self.para_xi.kernel.info.ks = 0;
            self.para_xi.kernel.info.kh = 1;
            self.para_xi.kernel.info.idx_sigma = [];
            self.para_xi.kernel.info.no_zero_peak = true;
            self.para_xi.kernel.info.df = self.df;
            self.para_xi.kernel.info.taperinfo = self.taperinfo;
        end
    otherwise
        error('no such method')
end

end
