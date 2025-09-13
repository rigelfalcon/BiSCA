function bc = calc_bc(self, bs, denominator)

    switch self.normalization
        case 'skewness'
            bc = sqrt(abs(bs) .^ 2 ./ (denominator .^ 2 + self.dn_reg)); %+self.dn_reg
        case {'haubrich'}
            bc = bs ./ ((denominator) + self.dn_reg); %+self.dn_reg

        case {'hagihira', 'hag'}
            bc = bs ./ ((denominator) + self.dn_reg);
        otherwise
            error('no such normalization')
    end

end