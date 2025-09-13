function NonlinearOptimization(self)
    freemu = isempty(self.mufix);

    if freemu
        self.mufix = self.para_alpha.kernel.para.mu;
    end

    fit_s(self);

    if freemu
        self.mufix = [];
    end
    fit_s(self);

    % get the theoretical spectrum and bispectrum with the lost function Brillinger 1985 and Leonenko proposed in 1998
    %%
    %%{
    if ~isempty(self.bs) && size(self.bs, 1) == size(self.s, 1) && size(self.bs, 2) == size(self.s, 1)
        self.para_bs.input_config = self.para_bs; %false

        self.para_bs.use_bssigma = false; %false
        self.para_bs.use_bsnu = false; %false
        self.para_bs.use_bsd = false; %false

        if self.para_bs.fix_mu
            self.mufix = self.para_alpha.kernel.para.mu;
        end

        fit_bs(self);

        % tune bispectrum
        self.para_bs.use_bssigma = self.para_bs.input_config.use_bssigma;
        self.para_bs.use_bsnu = self.para_bs.input_config.use_bsnu;
        self.para_bs.use_bsd = self.para_bs.input_config.use_bsd; %false

        if self.para_bs.use_bssigma || self.para_bs.use_bsnu || self.para_bs.use_bsd
            fit_bs(self);
        end

    else
        % spectra-only mode (no bispectrum data provided)
    end

    %%}
    %%
end