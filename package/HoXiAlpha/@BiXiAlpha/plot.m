function plot(self)

    sest = self.shat;
    subplot(2, 1, 1)
    plot(self.f, self.s);
    hold on

    subplot(2, 1, 2)
    plot(self.f, sum(sest, 2));
    hold on
end