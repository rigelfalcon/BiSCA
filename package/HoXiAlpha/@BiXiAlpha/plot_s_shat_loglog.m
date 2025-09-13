function plot_s_shat_loglog(self)

    pause(0.00001);
    hFig = gcf;

    FontSize = 6;
    LineWidth = 1;
    MarkerSize = 2;
    set(0, 'defaultAxesFontName', 'arial')
    set(0, 'defaultTextFontName', 'arial')

    p = get(0, "MonitorPositions");
    if size(p, 1) > 1
        Position = p(2, :);
    else
        Position = p(1, :);
    end

    Position = [Position(1) + 200, Position(2) + 200, 400, 150];
    hFig.Position = Position;

    %%

    cm = (plt.Set1);
    hold on;
    plot(NaN, 'DisplayName', 'real data', 'Color', cm(1, :));
    plot(log10(self.f), 10 * log10(mean(self.s, 2)), 'Color', cm(2, :), 'LineWidth', LineWidth, 'MarkerSize', MarkerSize);
    plot(log10(self.f), 10 * log10(self.shat), 'Color', cm(3, :), 'LineWidth', LineWidth, 'MarkerSize', MarkerSize);


    [~, scomp] = pred_s(self.f, self.ks, self.kh, self.para_xi, self.para_alpha);
    idx_drop = false(1, size(scomp, 2));
    idx_drop(2) = true;
    scomp(:, idx_drop) = [];

    if strcmpi(self.para_xi.type, 'tstudent+peak')
        scomp(:, 1) = scomp(:, 1) + scomp(:, 2);
        scomp(:, 2) = [];
    end

    smin = min(10 * [log10(mean(self.s, 2)); log10(self.shat)], [], 'all');
    smax = max(10 * [mean(log10(self.s), 2); log10(self.shat)], [], 'all');
    scomp(10 * log10(scomp) < smin - 10) = nan;

    for icomp = 1:size(scomp, 2) % plot(self.f, 10 * log10(scomp),'-', 'LineWidth', 2); hold on
        plot(log10(self.f), 10 * log10(scomp(:, icomp)), '-', 'Color', cm(3 + icomp, :), 'LineWidth', LineWidth);
        hold on

    end

    mu = get_harmonic_mu_sigma(self.para_alpha.kernel.para.mu, 1, self.ks, self.kh, self.para_alpha.no_zero_peak);
    mu(idx_drop(2:end)) = [];
    xline(mu, ':k', 'Alpha', 0.6, 'LineWidth', LineWidth, 'HandleVisibility', 'off')

    ylim([smin - 2, smax + 2])
    xticks([0:5:max(self.f)])
    xlabel('Frequency (Hz)')

    ylabel('Power Spectrum density (db)')

    hA = gca;
    hA.XAxis.FontSize = FontSize;
    hA.YAxis.FontSize = FontSize;
    hA.YRuler.TickLabelGapOffset = 5;
    xtickangle(45);

    if self.ks > 0
        lagalpha = [append('1/', string(num2str((self.ks + 1:-1:2)')), ' harmonic'); "alpha peak"; append(string(num2str((2:self.kh)')), '/1 harmonic')];
    else
        lagalpha = ["alpha peak"; append(string(num2str((2:self.kh)')), '/1 harmonic')];
    end

    if ~self.para_alpha.no_zero_peak
        lagalpha = ['0 peak'; lagalpha]; %\delta f peak
    end

    lagempty = repmat('', 1, self.ks + self.kh);
    ws = warning('off');
    lagalpha = lagalpha(1:self.ks + self.kh)';
    lagalpha(idx_drop(2:end)) = [];
    lgd = legend(["", "empirical spectra", "fitted spectra", "\xi process", lagalpha, lagempty], 'FontSize', FontSize, 'Location', 'eastoutside');
    lgd.FontName = 'Arial';
    warning(ws);

    drawnow;
    pause(0.05); % this innocent line prevents the Matlab hang

end