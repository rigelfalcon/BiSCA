function visualize_bc_stat_max_median_symlog_log(bc, Ns, p_threshold, panel_label)
    % VISUALIZE_BC_STAT_MAX_MEDIAN_SYMLOG_LOG Visualize bicoherence statistics
    %
    % R1-24 Changes:
    %   - Increased font sizes for better readability
    %   - Added optional panel_label parameter for (a), (b), etc.
    %   - Changed chi-square line from blue to magenta
    %
    % Usage:
    %   visualize_bc_stat_max_median_symlog_log(bc, Ns, p_threshold)
    %   visualize_bc_stat_max_median_symlog_log(bc, Ns, p_threshold, '(a)')
    
    if nargin < 4
        panel_label = '';
    end
    
    % Font size settings (tuned for subplot layout)
    FontSize_axis = 14;
    FontSize_label = 15;
    FontSize_text = 17;
    FontSize_legend = 14;
    FontSize_panel = 20;
    
    % Compute statistics for 'max' and 'median'
    stat_max = calc_bc_stat(bc, 'max', Ns, p_threshold);
    stat_median = calc_bc_stat(bc, 'median', Ns, p_threshold);

    % Prepare data: compute scaled bicoherence and remove NaNs
    bcvec = abs(reshape(bc, [], size(bc, 3))) .^ 2;
    scaled_bcvec = 2 * Ns * bcvec;
    scaled_bcvec = scaled_bcvec(~isnan(scaled_bcvec) & scaled_bcvec > 0);

    % Create figure
    ax = gca;

    % Compute histogram in symlog scale
    min_val = max(min(scaled_bcvec), 1e-10);
    max_val = max(scaled_bcvec);
    C = 1; % Symlog scaling constant
    z_min = log10(1 + min_val / C);
    z_max = log10(1 + max_val / C);
    num_bins = 50;
    z_edges = linspace(z_min, z_max, num_bins + 1);
    x_edges = C * (10 .^ z_edges - 1);
    [counts, ~] = histcounts(scaled_bcvec, x_edges);
    bin_widths = diff(x_edges);
    total_counts = sum(counts);
    pdf_heights = counts ./ (total_counts .* bin_widths);
    pdf_heights = max(pdf_heights, 1e-10); % Ensure positive for log scale

    % Plot histogram as a single patch
    hold on;
    base_value = 1e-10;
    x_vertices = [];
    y_vertices = [];

    for i = 1:length(pdf_heights)
        x_vertices = [x_vertices, x_edges(i), x_edges(i), x_edges(i + 1), x_edges(i + 1)];
        y_vertices = [y_vertices, base_value, pdf_heights(i), pdf_heights(i), base_value];
    end

    patch(x_vertices, y_vertices, [0.6 0.8 1], 'FaceAlpha', 0.6, 'EdgeColor', 'none', ...
        'DisplayName', 'Bicoherence Histogram');

    % Chi-square(2) PDF
    x_orig = linspace(min_val, max_val, 1000);
    pdf_chi2 = chi2pdf(x_orig, 2);
    pdf_chi2 = max(pdf_chi2, 1e-10);
    plot(x_orig, pdf_chi2, '-', 'LineWidth', 3, 'Color', [0.8, 0, 0.8], 'DisplayName', 'Chi-square(2)');  % R1-24: Changed from blue (hard to see) to magenta

    % Non-central chi-square PDF
    sample_mean = median(scaled_bcvec);
    lambda0 = max(sample_mean - 2, 0);
    pdf_ncx2 = ncx2pdf(x_orig, 2, lambda0);
    pdf_ncx2 = max(pdf_ncx2, 1e-10);
    plot(x_orig, pdf_ncx2, '-', 'LineWidth', 3, 'Color', [1 0.5 0], ...
        'DisplayName', sprintf('Non-central Chi-square\n(2, \\lambda_0 = %.2f)', lambda0));

    % Set left y-axis
    ylim([1e-4, 1]);
    set(ax, 'YScale', 'log');
    ylabel('Density (Bicoherence, log scale)', 'FontSize', FontSize_label);

    % Max statistic PDF
    x_max_orig = linspace(stat_max.pd.icdf(1e-10), stat_max.pd.icdf(1 - 1e-8), 1000);
    pdf_max = stat_max.pd.pdf(x_max_orig);
    pdf_max = max(pdf_max, 1e-10);
    plot(x_max_orig, pdf_max, '-', 'LineWidth', 3, 'DisplayName', 'Max Statistic PDF', 'Color', [0 0 0 0.3]);
    hold on;

    % Median statistic PDF
    x_median_orig = linspace(stat_median.pd.icdf(1e-12), stat_median.pd.icdf(1 - 1e-12), 1000);
    pdf_median = stat_median.pd.pdf(x_median_orig);
    pdf_median = max(pdf_median, 1e-10);
    plot(x_median_orig, pdf_median, 'Color', [0.5 0 0 0.3], 'LineStyle', '-', 'LineWidth', 3, ...
        'DisplayName', 'Median Statistic PDF');

    ylim_right = get(ax, 'YLim');
    % Vertical lines with corrected legend names
    % R1-24 fix: Stagger text labels at different heights to avoid overlap
    % Use 90%, 80%, 70%, 60% of plot height for the 4 labels
    % Heights in log-space fractions (y-axis is log scale, so linear fractions land near top)
    yl = ylim_right;
    log_min = log10(yl(1)); log_max = log10(yl(2));
    log_fracs = [0.70, 0.58, 0.46, 0.34];  % visual fractions in log space
    text_heights_abs = 10 .^ (log_min + log_fracs * (log_max - log_min));
    
    if stat_max.T > stat_max.th
        pos_x_max_stat = apply_symlog_offset(stat_max.T, 1, 1);
        pos_x_max_th = apply_symlog_offset(stat_max.th, -1, 1);
    else
        pos_x_max_stat = apply_symlog_offset(stat_max.T, -1, 1);
        pos_x_max_th = apply_symlog_offset(stat_max.th, 1, 1);
    end

    plot([stat_max.T, stat_max.T], [ylim_right(1), ylim_right(2)], 'r-.', 'LineWidth', 3, ...
        'DisplayName', 'Max Statistic');
    text(pos_x_max_stat, text_heights_abs(1), ...
        sprintf('%.3f', stat_max.bc_T), 'Rotation', 90, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', FontSize_text - 2, 'Color', 'r');

    plot([stat_max.th, stat_max.th], [ylim_right(1), ylim_right(2)], 'm--', 'LineWidth', 3, ...
        'DisplayName', 'Max Threshold');
    text(pos_x_max_th, text_heights_abs(2), ...
        sprintf('%.3f', stat_max.bc_th), 'Rotation', 90, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', FontSize_text - 2, 'Color', 'm');

    if stat_median.T > stat_median.th
        pos_x_median_stat = apply_symlog_offset(stat_median.T, 1, 1);
        pos_x_median_th = apply_symlog_offset(stat_median.th, -1, 1);
    else
        pos_x_median_stat = apply_symlog_offset(stat_median.T, -1, 1);
        pos_x_median_th = apply_symlog_offset(stat_median.th, 1, 1);
    end

    plot([stat_median.T, stat_median.T], [ylim_right(1), ylim_right(2)], '-.', 'LineWidth', 3, ...
        'Color', [0 0.6 0], 'DisplayName', 'Median Statistic');
    text(pos_x_median_stat, text_heights_abs(3), ...
        sprintf('%.3f', stat_median.bc_T), 'Rotation', 90, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', FontSize_text - 2, 'Color', [0 0.6 0]);

    plot([stat_median.th, stat_median.th], [ylim_right(1), ylim_right(2)], 'b--', 'LineWidth', 3, ...
        'DisplayName', 'Median Threshold');
    text(pos_x_median_th, text_heights_abs(4), ...
        sprintf('%.3f', stat_median.bc_th), 'Rotation', 90, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', FontSize_text - 2, 'Color', 'b');

    ylabel('Density (Statistic, log scale)', 'FontSize', FontSize_label);

    % Apply symlog transformation to x-axis
    symlogax(ax, 'x');
    xlims = xlim;
    xlim([xlims(1), xlims(2) + 0.05 * (xlims(2) - xlims(1))])

    % Customize plot (legend handled externally as shared legend)
    xlabel('t-statistic of Bicoherence (symlog scale)', 'FontSize', FontSize_label);
    grid on;
    set(ax, 'FontName', 'Helvetica', 'FontSize', FontSize_axis, 'GridAlpha', 0.3);

    % Add significance text - R1-24 fix: Use more compact format to avoid overlap
    if stat_max.h
        sig_text_max = 'Sig';
    else
        sig_text_max = 'NS';
    end

    if stat_median.h
        sig_text_median = 'Sig';
    else
        sig_text_median = 'NS';
    end

    % Place text in bottom-left corner with compact format
    text(0.02, 0.02, sprintf('Max p=%.4f (%s) | Med p=%.4f (%s) | \\alpha=%.4f', ...
        stat_max.p, sig_text_max, stat_median.p, sig_text_median, p_threshold), ...
        'Units', 'normalized', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
        'FontSize', FontSize_legend - 1, 'BackgroundColor', [1 1 1 0.8]);
    
    % R1-24: Add panel label (a), (b), etc. if provided
    % Position adjusted to avoid overlap with Y-axis label
    if ~isempty(panel_label)
        text(0.15, 0.98, panel_label, 'Units', 'normalized', ...
            'FontSize', FontSize_panel, 'FontWeight', 'bold', ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
            'BackgroundColor', [1 1 1 0.7]);
    end
end

function symlogax(ax, axis)
    % Apply symmetric logarithmic scale to the specified axis
    C = 1; % Scaling constant

    % Transform objects with XData (Lines and Patches)
    h = findobj(ax, '-property', 'XData');

    for i = 1:length(h)

        switch lower(axis)
            case 'x'
                x = get(h(i), 'XData');
                set(h(i), 'XData', sign(x) .* log10(1 + abs(x) / C));
            case 'y'
                y = get(h(i), 'YData');
                set(h(i), 'YData', sign(y) .* log10(1 + abs(y) / C));
        end

    end

    % Transform text objects' x-positions
    text_objs = findobj(ax, 'Type', 'Text');

    for i = 1:length(text_objs)
        pos = get(text_objs(i), 'Position');

        switch lower(axis)
            case 'x'
                x = pos(1); % x-coordinate of text position
                new_x = sign(x) .* log10(1 + abs(x) / C);
                set(text_objs(i), 'Position', [new_x, pos(2), pos(3)]);
            case 'y'
                y = pos(2); % y-coordinate of text position
                new_y = sign(y) .* log10(1 + abs(y) / C);
                set(text_objs(i), 'Position', [pos(1), new_y, pos(3)]);
        end

    end

    % Set custom tick marks
    switch lower(axis)
        case 'x'
            xlim = get(ax, 'XLim');
            t0 = ceil(log10(C)):ceil(log10(max(abs(xlim))));
            t1 = 10 .^ t0;
            ticks = [0, t1];
            set(ax, 'XTick', symlog_transform(ticks, C));
            set(ax, 'XTickLabel', arrayfun(@(x) sprintf('%.0f', x), ticks, 'UniformOutput', false));
        case 'y'
            ylim = get(ax, 'YLim');
            t0 = ceil(log10(C)):ceil(log10(max(abs(ylim))));
            t1 = 10 .^ t0;
            ticks = [0, t1];
            set(ax, 'YTick', symlog_transform(ticks, C));
            set(ax, 'YTickLabel', arrayfun(@(x) sprintf('%.0f', x), ticks, 'UniformOutput', false));
    end

end

function y = symlog_transform(x, C)
    y = sign(x) .* log10(1 + abs(x) / C);
end

function x_corrected = apply_symlog_offset(x, offset_dir, C)
    % Apply an offset in symlog space and return the corrected position in data space
    % Inputs:
    %   x: Position in data space (e.g., stat_max.T)
    %   offset_dir: +1 for right offset, -1 for left offset in symlog space
    %   C: Scaling constant for symlog transformation (default=1)
    % Output:
    %   x_corrected: Corrected position in data space
    if offset_dir < 1
        offset_symlog = 0.01; % Fixed offset in symlog space
    else
        offset_symlog = 0.12; % Fixed offset in symlog space
    end

    % Transform to symlog space
    x_symlog = sign(x) .* log10(1 + abs(x) / C);
    % Apply offset
    x_symlog_corrected = x_symlog + offset_dir * offset_symlog;
    % Back-transform to data space
    x_corrected = sign(x_symlog_corrected) .* (10 .^ abs(x_symlog_corrected) - 1) * C;
end