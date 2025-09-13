function visualize_bc_stat_median(bc, Ns, p_threshold)
    % Compute statistics for 'median'
    stat_median = calc_bc_stat(bc, 'median', Ns, p_threshold);

    % Prepare data: compute scaled bicoherence and remove NaNs
    bcvec = abs(reshape(bc, [], size(bc, 3))) .^ 2;
    scaled_bcvec = 2 * Ns * bcvec;
    scaled_bcvec = scaled_bcvec(~isnan(scaled_bcvec) & scaled_bcvec > 0);

    % Create figure
    ax = gca;

    % Compute histogram bins and heights manually for symlog transformation
    min_val = max(min(scaled_bcvec), 1e-10);
    max_val = max(scaled_bcvec);
    edges = linspace(min_val, max_val, 50); % Linear bins in original space
    [counts, bin_edges] = histcounts(scaled_bcvec, edges, "Normalization", "pdf");
    bin_width = diff(bin_edges);
    pdf_heights = counts / (sum(counts) * bin_width(1)); % Normalize to PDF
    pdf_heights = max(pdf_heights, 1e-10); % Ensure positive for log scale

    % Plot histogram as a single patch object with multiple rectangles
    hold on;
    base_value = 1e-10; % Small positive value for log scale
    x_vertices = [];
    y_vertices = [];
    for i = 1:length(pdf_heights)
        % Define vertices for each bin
        x_vertices = [x_vertices, bin_edges(i), bin_edges(i), bin_edges(i + 1), bin_edges(i + 1)];
        y_vertices = [y_vertices, base_value, pdf_heights(i), pdf_heights(i), base_value];
    end
    patch(x_vertices, y_vertices, [0.6 0.8 1], 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Data Histogram');

    % Chi-square(2) PDF
    x_orig = linspace(min_val, max_val, 1000);
    pdf_chi2 = chi2pdf(x_orig, 2);
    pdf_chi2 = max(pdf_chi2, 1e-10); % Ensure positive for log scale
    plot(x_orig, pdf_chi2, '-', 'LineWidth', 3, 'Color', [0, 0, 1, 0.2], 'DisplayName', 'Chi-square(2)');

    % Compute and plot non-central chi-square PDF
    sample_mean = median(scaled_bcvec);
    lambda0 = max(sample_mean - 2, 0); % Estimate lambda0 as in original code
    pdf_ncx2 = ncx2pdf(x_orig, 2, lambda0);
    pdf_ncx2 = max(pdf_ncx2, 1e-10); % Ensure positive for log scale

    % Set left y-axis to log scale
    ylim([1e-4, 1]);
    set(ax, 'YScale', 'log');
    ylabel('Density (Bicoherence, log scale)', 'FontSize', 14);

    % Median statistic PDF
    x_median_orig = linspace(stat_median.pd.icdf(1e-12), stat_median.pd.icdf(1 - 1e-12), 1000);
    pdf_median = stat_median.pd.pdf(x_median_orig);
    pdf_median = max(pdf_median, 1e-10); % Ensure positive for log scale
    plot(x_median_orig, pdf_median, 'Color', [0.5 0 0 0.3], 'LineStyle', '-', 'LineWidth', 3, 'DisplayName', ['PDF of median bicoherence']);
    hold on;

    % Vertical lines for median statistic and threshold with vertical text
    ylim_right = get(ax, 'YLim');
    if stat_median.T > stat_median.th
        pos_x_median_stat = apply_symlog_offset(stat_median.T, 1, 1);
        pos_x_median_th = apply_symlog_offset(stat_median.th, -1, 1);
    else
        pos_x_median_stat = apply_symlog_offset(stat_median.T, -1, 1);
        pos_x_median_th = apply_symlog_offset(stat_median.th, 1, 1);
    end
    plot([stat_median.T, stat_median.T], [ylim_right(1), ylim_right(2)], 'g-.', 'LineWidth', 3, 'DisplayName', 'Sample test statistic');
    text(pos_x_median_stat, ylim_right(2), 'Sample test statistic', 'Rotation', 90, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 14);
    plot([stat_median.th, stat_median.th], [ylim_right(1), ylim_right(2)], 'b--', 'LineWidth', 3, 'DisplayName', 'Critical value');
    text(pos_x_median_th, ylim_right(2), 'Critical value', 'Rotation', 90, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'FontSize', 14);

    ylabel('Density (Statistic, log scale)', 'FontSize', 12);

    % Apply symlogax to x-axis
    symlogax(ax, 'x');

    % Customize plot aesthetics
    xlabel('t-statistic of Bicoherence (symlog scale)', 'FontSize', 12);
    legend('Location', 'eastoutside', 'FontSize', 10);
    grid on;
    set(ax, 'FontName', 'Helvetica', 'FontSize', 14, 'GridAlpha', 0.3);

    xlims = xlim;
    xlim([xlims(1), xlims(2) + 0.1])

    % Add significance text
    if stat_median.h
        sig_text_median = 'Significant';
    else
        sig_text_median = 'Not Significant';
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
    if offset_dir < 0
        offset_symlog = 0.01; % Fixed offset in symlog space
    else
        offset_symlog = 0.12; % Fixed offset in symlog space
    end
    % Transform to symlog space
    x_symlog = sign(x) .* log10(1 + abs(x) / C);
    % Apply offset
    x_symlog_corrected = x_symlog + offset_dir * offset_symlog;
    % Back-transform to data space
    x_corrected = sign(x_symlog_corrected) .* (10.^abs(x_symlog_corrected) - 1) * C;
end