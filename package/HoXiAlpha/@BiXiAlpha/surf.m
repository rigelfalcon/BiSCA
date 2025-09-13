function surf(self)
    bhat = pred_bs(self);
    bsre = real(bhat);
    bsim = imag(bhat);

    subplot(2, 1, 1)
    p1 = surf(self.fx, self.fy, abs(self.bs)); hold on;
    view(0, 90)
    set(p1, 'edgecolor', 'none')
    shading interp
    grid off
    colormap(flip(hot))
    colorbar

    subplot(2, 1, 2)
    p2 = surf(self.fx, self.fy, sqrt(bsre .^ 2 + bsim .^ 2));
    view(0, 90)
    set(p2, 'edgecolor', 'none')
    shading interp
    grid off
    colormap(flip(hot))
    colorbar
end
