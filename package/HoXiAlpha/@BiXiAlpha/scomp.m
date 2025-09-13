function s_comp = scomp(self)
[~, s_comp] = pred_s(self.f, self.ks, self.kh, self.para_xi, self.para_alpha);
end
