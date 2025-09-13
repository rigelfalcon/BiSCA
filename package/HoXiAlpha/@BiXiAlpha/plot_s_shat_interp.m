function plot_s_shat_interp(self)

    pause(0.00001);
    hFig = gcf;

    FontSize = 6;
    LineWidth = 2;
    MarkerSize = 2;
    set(0, 'defaultAxesFontName', 'arial')
    set(0, 'defaultTextFontName', 'arial')

    p = get(0, "MonitorPositions");
    if size(p, 1) > 1
        Position = p(2, :);
    else
        Position = p(1, :);
    end

    Position = [Position(1) + 200, Position(2) + 200, 300, 180];
    hFig.Position = Position;

    %%
    f = self.f;
    s = self.s;
    shat = self.shat;

    h = 0.3; %0.4

    f_interp = linspace(min(f), max(f), length(f) * 10)';
    [s] = NwSmoothInline(f, s, h, f_interp);
    [shat] = NwSmoothInline(f, shat, h, f_interp);

    [~, scomp] = pred_s(f_interp, self.ks, self.kh, self.para_xi, self.para_alpha);
    f = f_interp;

    %%

    cm = (plt.Set1);
    hold on;
    plot(NaN, 'DisplayName', 'real data', 'Color', cm(1, :));
    plot(f, 10 * log10(mean(s, 2)), 'Color', cm(2, :), 'LineWidth', LineWidth, 'MarkerSize', MarkerSize);
    plot(f, 10 * log10(shat), 'Color', cm(3, :), 'LineWidth', LineWidth, 'MarkerSize', MarkerSize);

    idx_drop = false(1, size(scomp, 2));
    idx_drop(2) = true;
    scomp(:, idx_drop) = [];

    if strcmpi(self.para_xi.type, 'tstudent+peak')
        scomp(:, 1) = scomp(:, 1) + scomp(:, 2);
        scomp(:, 2) = [];
    end

    smin = min(10 * [log10(mean(s, 2)); log10(shat)], [], 'all');
    smax = max(10 * [mean(log10(s), 2); log10(shat)], [], 'all');
    scomp(10 * log10(scomp) < smin - 10) = nan;

    for icomp = 1:size(scomp, 2) % plot(self.f, 10 * log10(scomp),'-', 'LineWidth', 2); hold on
        plot(f, 10 * log10(scomp(:, icomp)), '-', 'Color', cm(3 + icomp, :), 'LineWidth', LineWidth);
        hold on

    end

    if strcmpi(self.peak_relation, 'harmonic')
        mu = get_harmonic_mu_sigma(self.para_alpha.kernel.para.mu, 1, self.ks, self.kh, self.para_alpha.no_zero_peak);
    elseif strcmpi(self.peak_relation, 'free')
        mu = self.para_alpha.kernel.para.mu(:)';
    end
    xline(mu, ':k', 'Alpha', 0.6, 'LineWidth', LineWidth / 2, 'HandleVisibility', 'off')

    xlim([0, max(f)])
    ylim([smin - 2, smax + 2])
    xticks([0:5:max(f)])
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
    warning(ws);

    drawnow;
    pause(0.05); % this innocent line prevents the Matlab hang

end

function [yq, L, dbg] = NwSmoothInline(x, y, h, xq)
    if nargin < 4 || isempty(xq)
        xq = x;
    end
    [n, dx] = size(x);
    [nq, ~] = size(xq);
    if size(x, 2) > 1
        x = reshape(x, [n, 1, dx]);
    end
    if size(xq, 2) > 1
        xq = reshape(xq, [nq, 1, dx]);
    end
    D = x - permute(xq, [2, 1, 3]);
    Dkn = gaussian_kernel(D, h);
    yq = (sum(Dkn .* y, 1) ./ sum(Dkn, 1))';
    dbg.s = sum(Dkn, 1);
    if nargout > 1
        L = Dkn ./ sum(Dkn, 1).';
    end
end

function u = gaussian_kernel(u, b)
    u = u ./ b;
    u = (1 / sqrt(2 * pi)) * exp(-0.5 * sum((u .^ 2), 3));
end