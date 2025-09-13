function shat=get_shat_component(self,component)
    switch component
        case {'xi'}
            shat=pred_s_xi(self.f,self.para_xi);
        case {'alpha'}
            shat=pred_s_alpha(self.f,self.kh,self.para_alpha);
        otherwise
            error('the component is not exsit, only xi and alpha');
    end

end