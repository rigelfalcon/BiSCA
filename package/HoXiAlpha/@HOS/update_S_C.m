function update_S_C(self)
    if self.update_cross && ~self.update_lock
        if strcmpi(self.method, 'fft')
            U2 = size(self.coef, 2);
        elseif strcmpi(self.method, 'pmtm')
            U2 = size(self.coef, 2) ./ self.window;
        else
            error('not support')
        end
        X = permute(self.coef, [3, 2, 1]); %#ok<*PROP>
        self.Sloc = pagemtimes(X, "none", X, "ctranspose") ./ U2;
        ss = sqrt(permute(self.s, [2, 3, 1]));
        ss = pagemtimes(ss, 'none', ss, 'ctranspose'); %toc% tic;ss=ss.*permute(ss,[2,1,3]);toc% ss=diag2full(self.s');
        self.Cloc = self.Sloc ./ ss; %[coh,psd]=calc_cov2corr(self.Sloc)
    end

end