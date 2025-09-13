function plot_s_shat_spa(self)

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

    Position = [200 200 293.3333 149.3333];
    hFig.Position = Position;

    %%

    % IBM colorblind-safe palette 
    % Peaks use colors that avoid conflict with inferno/viridis colormaps.
    % Non-peak lines use warm colors (Orange/Yellow/Vermillion).
    
    % Lines NOT drawn on bicoherence plots (can use warm colors)
    c_data  = [0.90, 0.60, 0.00];  % Multitaper spectrum -> Orange
    c_fit   = [0.00, 0.45, 0.70];  % SPA spectrum -> IBM Blue
    c_xi    = [0.00, 0.62, 0.45];  % xi process -> Bluish Green (IBM)
    
    % Peak lines drawn on bicoherence plots (avoid inferno/viridis overlap)
    c_alpha = [1.00, 0.00, 0.00];  % alpha peak -> Pure Red
    c_h2    = [0.60, 0.00, 1.00];  % 2/1 harmonic -> Bright Violet
    c_h3    = [1.00, 1.00, 0.00];  % 3/1 harmonic -> Bright Yellow

    hold on;
    h_data = plot(self.f, 10 * log10(mean(self.s, 2)), 'Color', c_data, 'LineWidth', LineWidth, 'MarkerSize', MarkerSize, 'DisplayName', 'Multitaper spectrum');
    h_fit  = plot(self.f, 10 * log10(self.shat), 'Color', c_fit,  'LineWidth', LineWidth, 'MarkerSize', MarkerSize, 'DisplayName', 'SPA spectrum');

    [~, scomp] = pred_s(self.f, self.ks, self.kh, self.para_xi, self.para_alpha);

    if strcmpi(self.para_xi.type, 'tstudent+peak')
        scomp(:, 1) = scomp(:, 1) + scomp(:, 2);
        scomp(:, 2) = []; % This modifies scomp and its number of columns
    end

    smin = min(10 * [log10(mean(self.s, 2)); log10(self.shat)], [], 'all');
    smax = max(10 * [mean(log10(self.s), 2); log10(self.shat)], [], 'all'); % Corrected log10 placement for smax
    scomp(10 * log10(scomp) < smin - 10) = nan;

    h_comp = gobjects(1, size(scomp, 2));
    for icomp = 1:size(scomp, 2)
        % In this plot, scomp(:,1) is expected to be the \xi component; peaks follow.
        if icomp == 1
            cl = c_xi;
            dname = "\\xi process";
        elseif icomp == 2
            cl = c_alpha;
            dname = "alpha peak";
        elseif icomp == 3
            cl = c_h2;
            dname = "2/1 harmonic";
        elseif icomp == 4
            cl = c_h3;
            dname = "3/1 harmonic";
        else
            cl = c_fit; % de-emphasize extra components if present
            dname = "peak " + string(icomp - 1);
        end
        h_comp(icomp) = plot(self.f, 10 * log10(scomp(:, icomp)), '-', 'Color', cl, 'LineWidth', LineWidth, 'DisplayName', dname);
        hold on
    end

    mu = get_harmonic_mu_sigma(self.para_alpha.kernel.para.mu(self.ks + self.k0 + 1), 1, self.ks, self.kh, self.para_alpha.no_zero_peak);
    xline(mu, ':k', 'Alpha', 0.6, 'LineWidth', LineWidth, 'HandleVisibility', 'off')

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

    ws = warning('off');
    lgd = legend([h_data, h_fit, h_comp(:)'], 'FontSize', FontSize, 'Location', [0.553521337313275, 0.498243415474555, 0.380681818181818, 0.4453125]);
    lgd.FontName = 'Arial';
    warning(ws);

    drawnow;
    pause(0.05); % this innocent line prevents the Matlab hang

end