function plot_s_shat_struct(bxa_struct)
% Wrapper for plot_s_shat that accepts a precomputed struct instead of BiXiAlpha object
% Usage: plot_s_shat_struct(bixialpha)  where bixialpha is a struct with precomputed fields

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

    % IBM colorblind-safe palette 
    c_data  = [0.90, 0.60, 0.00];  % Multitaper spectrum -> Orange
    c_fit   = [0.00, 0.45, 0.70];  % BiSCA spectrum -> IBM Blue
    c_xi    = [0.00, 0.62, 0.45];  % xi process -> Bluish Green (IBM)
    c_alpha = [1.00, 0.00, 0.00];  % alpha peak -> Pure Red
    c_h2    = [0.60, 0.00, 1.00];  % 2/1 harmonic -> Bright Violet
    c_h3    = [1.00, 1.00, 0.00];  % 3/1 harmonic -> Bright Yellow

    % Check required fields
    required = {'f', 's', 'shat', 'ks', 'kh', 'k0', 'para_xi', 'para_alpha'};
    for i = 1:numel(required)
        if ~isfield(bxa_struct, required{i})
            error('Missing required field: %s', required{i});
        end
    end

    hold on;
    h_data = plot(bxa_struct.f, 10 * log10(mean(bxa_struct.s, 2)), 'Color', c_data, 'LineWidth', LineWidth, 'MarkerSize', MarkerSize, 'DisplayName', 'Multitaper spectrum');
    h_fit  = plot(bxa_struct.f, 10 * log10(bxa_struct.shat), 'Color', c_fit,  'LineWidth', LineWidth, 'MarkerSize', MarkerSize, 'DisplayName', 'BiSCA spectrum');

    [~, scomp] = pred_s(bxa_struct.f, bxa_struct.ks, bxa_struct.kh, bxa_struct.para_xi, bxa_struct.para_alpha);

    if strcmpi(bxa_struct.para_xi.type, 'tstudent+peak')
        scomp(:, 1) = scomp(:, 1) + scomp(:, 2);
        scomp(:, 2) = [];
    end

    smin = min(10 * [log10(mean(bxa_struct.s, 2)); log10(bxa_struct.shat)], [], 'all');
    smax = max(10 * [log10(mean(bxa_struct.s, 2)); log10(bxa_struct.shat)], [], 'all');
    scomp(10 * log10(scomp) < smin - 5) = nan;
    idx_drop = sum(scomp, 1, "omitmissing") == 0;
    scomp(:, idx_drop) = [];

    h_comp = gobjects(1, size(scomp, 2));
    for icomp = 1:size(scomp, 2)
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
            cl = c_fit;
            dname = "peak " + string(icomp - 1);
        end
        h_comp(icomp) = plot(bxa_struct.f, 10 * log10(scomp(:, icomp)), '-', 'Color', cl, 'LineWidth', LineWidth, 'DisplayName', dname);
        hold on
    end

    if length(bxa_struct.para_alpha.kernel.para.mu) == bxa_struct.ks + bxa_struct.k0 + bxa_struct.kh
        mu = get_harmonic_mu_sigma(bxa_struct.para_alpha.kernel.para.mu(bxa_struct.ks + bxa_struct.k0 + 1), 1, bxa_struct.ks, bxa_struct.kh, bxa_struct.para_alpha.no_zero_peak);
    else
        mu = get_harmonic_mu_sigma(bxa_struct.para_alpha.kernel.para.mu, 1, bxa_struct.ks, bxa_struct.kh, bxa_struct.para_alpha.no_zero_peak);
    end
    mu(idx_drop(2:end)) = [];
    xline(mu, ':k', 'Alpha', 0.6, 'LineWidth', LineWidth, 'HandleVisibility', 'off')

    xlim([0, max(bxa_struct.f)])
    ylim([smin - 1, smax + 1])
    xticks([0:5:max(bxa_struct.f)])
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
    pause(0.05);

end
