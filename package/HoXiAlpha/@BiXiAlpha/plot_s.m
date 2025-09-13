function plot_s_shat(self)

    pause(0.00001);
    hFig = gcf;

    FontSize = 6;
    LineWidth = 1;
    MarkerSize = 2;
    set(0, 'defaultAxesFontName', 'arial')
    set(0, 'defaultTextFontName', 'arial')

    p = get(0, "MonitorPositions");
    Position = p(1, :); %p(2,:);
    Position = [Position(1) + 200, Position(2) + 200, 270, 150];
    hFig.Position = Position;

    %%

    cm = plt.tab10;
    hold on;
    plot(NaN, 'DisplayName', 'real data', 'Color', cm(1, :));
    plot(self.f, 10 * log10(mean(self.s, 2)), 'Color', cm(2, :), 'LineWidth', LineWidth, 'MarkerSize', MarkerSize);

    smin = min(10 * [log10(mean(self.s, 2)); log10(self.shat)], [], 'all');
    smax = max(10 * [mean(log10(self.s), 2); log10(self.shat)], [], 'all');

    xlim([0, max(self.f)])
    ylim([smin - 2, smax + 2])
    xticks([0:5:max(self.f)])
    xlabel('Frequency (Hz)')

    ylabel('Spectrum (dB re 1 a.u.$^2$/Hz)', 'Interpreter', 'latex')

    hA = gca;
    hA.XAxis.FontSize = FontSize;
    hA.YAxis.FontSize = FontSize;
    hA.YRuler.TickLabelGapOffset = 5;
    xtickangle(45);

    lgd = legend(["", "empirical spectra"], 'FontSize', FontSize);
    lgd.FontName = 'Arial';

    drawnow;
    pause(0.05); % this innocent line prevents the Matlab hang

end