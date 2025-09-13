function plot_bsgrid(range)
    % range: The maximum value for the x:axis (e.g., max(self.f))

    hold on;
    c = [1, 0, 0];
    range_f1_f2 = [range(:); range(:) + max(range(:)) + range(2) - range(1)]';
    xline(range, 'Color', c, 'LineStyle', ':', 'LineWidth', 0.2);
    yline(range, 'Color', c, 'LineStyle', ':', 'LineWidth', 0.2);

    z = (zlim);
    z = z(2);

    for i = range_f1_f2
        h = plot3([0, i], [i, 0], [z, z], 'Color', c, 'LineStyle', ':', 'LineWidth', 0.2); %,'Layer', 'top' zorder
        uistack(h, 'top');
    end

end