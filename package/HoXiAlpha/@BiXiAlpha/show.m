function show(self)

    pause(0.00001);
    hFig = gcf;
    clf
    hFig.WindowState = 'maximized';

    if ~isempty(self.bs)
        subplot(3, 1, 1)
    end

    cm = colormap('lines');
    hold on;
    plot(NaN, 'DisplayName', 'real data', 'Color', cm(1, :));
    plot(self.f, 10 * mean(log10(self.s), 2), '-x', 'Color', cm(2, :), 'LineWidth', 2);
    plot(self.f, 10 * log10(self.shat), '-o', 'Color', cm(3, :), 'LineWidth', 2);

    [~, scomp] = pred_s(self.f, self.ks, self.kh, self.para_xi, self.para_alpha);

    if strcmpi(self.para_xi.type, 'tstudent+peak')
        scomp(:, 1) = scomp(:, 1) + scomp(:, 2);
        scomp(:, 2) = [];
    end

    smin = min(10 * [mean(log10(self.s), 2); log10(self.shat)], [], 'all');
    smax = max(10 * [mean(log10(self.s), 2); log10(self.shat)], [], 'all');

    for icomp = 1:size(scomp, 2)
        plot(self.f, 10 * log10(scomp(:, icomp)), '-', 'Color', cm(3 + icomp, :), 'LineWidth', 2);
        hold on
    end
    if strcmpi(self.peak_relation, 'harmonic')
        mu = get_harmonic_mu_sigma(self.para_alpha.kernel.para.mu, 1, self.ks, self.kh, self.para_alpha.no_zero_peak);
    elseif strcmpi(self.peak_relation, 'free')
        mu = self.para_alpha.kernel.para.mu(:)';
    end
    xline(mu, ':k', 'Alpha', 0.6, 'LineWidth', 5, 'HandleVisibility', 'off')

    xlim([0, max(self.f)])
    ylim([smin - 2, smax + 2])
    xticks([0:1:max(self.f)])

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
    legend(["real data", "mean of data", "fitted data", "\xi process", lagalpha(1:self.ks + self.kh)', lagempty])
    warning(ws);
    title('spectrum')

    if isempty(self.bs)
        return;
    end

    %% third order part
    bc = mean(self.bc, 3);
    bchat = mean(self.bchat, 3);
    bc_xi = mean(self.bchat_component('xi'), 3);
    bc_alpha = mean(self.bchat_component('alpha'), 3);

    dtick = 5;
    range = [0:dtick:max(self.f)];

    subplot(3, 4, 5)
    surfbc(self.fx, self.fy, abs(bc));
    xticks(range)
    yticks(range)
    hold on;
    plot_bsgrid(range)
    title(['orignal bicoherence [', self.normalization, ']'])

    subplot(3, 4, 9)
    surfbc(self.fx, self.fy, abs(bchat));
    xticks(range)
    yticks(range)
    hold on;
    plot_bsgrid(range)
    title(['fitted bicoherence [', self.normalization, ']'])

    subplot(3, 4, 6)
    surfbs(self.fx, self.fy, real(bc));
    xticks(range)
    yticks(range)
    hold on;
    plot_bsgrid(range)
    title('real part of orignal bicoherence')

    subplot(3, 4, 10)
    surfbs(self.fx, self.fy, real(bchat));
    xticks(range)
    yticks(range)
    hold on;
    plot_bsgrid(range)
    title('real part of fitted bicoherence')

    subplot(3, 4, 7)
    surfbs(self.fx, self.fy, imag(bc));
    xticks(range)
    yticks(range)
    hold on;
    plot_bsgrid(range)
    title('imaginary part of orignal bicoherence')

    subplot(3, 4, 11)
    surfbs(self.fx, self.fy, imag(bchat));
    xticks(range)
    yticks(range)
    hold on;
    plot_bsgrid(range)
    title('imaginary part of fitted bicoherence')

    %%
    %%{
    subplot(3, 4, 8)
    surfbs(self.fx, self.fy, abs(bc_xi));
    xticks(range)
    yticks(range)
    hold on;
    plot_bsgrid(range)
    title('bicoherence of \xi process')

    subplot(3, 4, 12)
    surfbc(self.fx, self.fy, abs(bc_alpha));
    xticks(range)
    yticks(range)
    hold on;
    plot_bsgrid(range)
    title('bicoherence of \alpha process')
    %%}

    if ~isempty(self.name)
        sgtitle(replace(append(self.name, ['-', self.hos_components]), '_', '-'))
    end

end