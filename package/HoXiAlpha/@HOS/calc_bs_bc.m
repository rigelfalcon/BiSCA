function calc_bs_bc(self)

    self.Xloc = self.Xloc - mean(self.Xloc, 1);
    switch lower(self.method)
        case {'hosa'}
            [self.bcloc, self.floc, self.bsloc, self.sloc] = call_bicoher(self.Xloc, self);
            self.wloc = self.floc / self.Fs * 2 * pi;
        case {'myhosa'}
            [self.bcloc, self.wloc, self.floc, self.bsloc, self.sloc] = calc_bicoher(self.Xloc, self.Fs, self.normalization, self.nfft, self.frange, self.compact, self.window, self.overlap);
        case {'fft'}

            [self.bcloc, self.wloc, self.floc, self.bsloc, self.sloc, self.bsdenomloc, self.coefloc, ~, self.bsvarloc, self.svarloc, self.sresloc, self.bsresloc] = calc_multiple_bicoher_noseg(self.Xloc, self.Fs, self.nfft, self.frange, self.compact, self.normalization, self.tf_return_seg, self.log_var_s);
            self.coefloc = reshape(self.coefloc, size(self.coefloc, 1), [], size(self.coefloc, 4));
            self.W2loc = [];
            self.W3loc = [];
            self.coefrawloc = [];
        case {'mtm'}
            [self.bcloc, self.wloc, self.floc, self.bsloc, self.sloc] = calc_bicoher_mtm(self.Xloc, self.Fs, self.nfft, self.frange, self.compact, false, [], self.NW, self.isseg, self.normalization);
        case {'pmtm'}

            if self.tf_return_res
                [self.bcloc, self.wloc, self.floc, self.bsloc, self.sloc, self.bsdenomloc, ~, ~, self.bsvarloc, self.svarloc, self.sresloc, self.bsresloc] = calc_multiple_bicoher_pmtm_noseg_large_sample_compact(self.Xloc, self.Fs, self.nfft, self.frange, self.compact, self.NW, self.normalization, self.taper, self.droplast, self.tf_return_seg, self.log_var_s); %,self.coefrawloc
            elseif self.tf_return_var
                [self.bcloc, self.wloc, self.floc, self.bsloc, self.sloc, self.bsdenomloc, ~, ~, self.bsvarloc, self.svarloc] = calc_multiple_bicoher_pmtm_noseg_large_sample_compact(self.Xloc, self.Fs, self.nfft, self.frange, self.compact, self.NW, self.normalization, self.taper, self.droplast, self.tf_return_seg, self.log_var_s); %,self.coefrawloc
            else
                [self.bcloc, self.wloc, self.floc, self.bsloc, self.sloc, self.bsdenomloc] = calc_multiple_bicoher_pmtm_noseg_large_sample_compact(self.Xloc, self.Fs, self.nfft, self.frange, self.compact, self.NW, self.normalization, self.taper, self.droplast, self.tf_return_seg, self.log_var_s); %,self.coefrawloc
            end

            if self.tf_return_seg
                % correct Ns to be the number of non-zero segments
                self.Ns = size(self.sloc, 2);
            end

            self.coefloc = reshape(self.coefloc, size(self.coefloc, 1), [], size(self.coefloc, 5));
            self.coefrawloc = [];
        case {'guido'}
            [self.bcloc, self.bsloc, self.floc, self.wloc, self.sloc] = call_compute_bispectrum(self.Xloc, self.Fs, self.nfft, self.frange, self.compact, self.window, self.overlap);

        otherwise
            error('no such method')
    end

end