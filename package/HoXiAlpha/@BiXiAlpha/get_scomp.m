function scomp=get_scomp(self)
[~, scomp] = pred_s(self.f, self.ks, self.kh, self.para_xi, self.para_alpha);
end