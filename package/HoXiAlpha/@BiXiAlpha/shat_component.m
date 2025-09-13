function shat = shat_component(self, component)

    switch component
        case {'xi'}
            [shat] = pred_s_xi(self.f, self.para_xi);
        case {'alpha'}
            [shat] = pred_s_alpha(self.f, self.ks, self.kh, self.para_alpha); %,comp gaussian_kernel(f,mu,sigma).*para_alpha.kernel.alpha.';

        otherwise
            error('the component is not exsit, only xi and alpha');
    end

end