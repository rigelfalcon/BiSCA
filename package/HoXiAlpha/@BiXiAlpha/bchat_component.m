function bchat = bchat_component(self, component, para_bs)

    if nargin < 3 || isempty(para_bs)
        para_bs = self.para_bs;
    end

    bshat = bshat_component(self, component, para_bs);

    switch component
        case {'xi', 'alpha'}
            bchat = calc_bc(self, bshat, self.bshatdenom);
        case {'xi_xi'}
            shat = pred_s_xi(self.f, self.para_xi);
            sxsysxyhat = get_sxsysxy(self.f, shat, [], self.fxy);
            bshatdenom = sqrt(sxsysxyhat);
            bchat = calc_bc(self, bshat, bshatdenom);
        case {'alpha_alpha'}
            shat = pred_s_alpha(self.f, self.ks, self.kh, self.para_alpha);
            sxsysxyhat = get_sxsysxy(self.f, shat, [], self.fxy);
            bshatdenom = sqrt(sxsysxyhat);
            bchat = calc_bc(self, bshat, bshatdenom);
        otherwise
            error('Invalid component for bchat computation.');
    end

end