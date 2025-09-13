function plot_bc_test(self, component, metric)

    hFig = gcf;
    FontSize = 16;
    LineWidth = 1;

    p = get(0, "MonitorPositions");
    if size(p, 1) > 1
        Position = p(2, :);
    else
        Position = p(1, :);
    end
    Position = [Position(1) + 200, Position(2) + 200, 550, 550];
    hFig.Position = Position;

    switch component
        case 'raw'
            bc = mean(self.bc, 3);
        case 'fit'
            bc = self.bchat;
        case 'res'
            bc = mean(self.bc, 3) - self.bchat; % "residual" (?)
        case 'xi'
            bc = self.bchat_component('xi');
        case 'alpha'
            bc = self.bchat_component('alpha');
        case 'xi_xi'
            bc = self.bchat_component('xi_xi');
        case 'alpha_alpha'
            bc = self.bchat_component('alpha_alpha');
        otherwise
            error('no such component')
    end

    T = abs(bc);
    Tmax = max(abs(mean(self.bc, 3)), [], "all");

    stat = self.stat.(component).(metric);

    stat.th = median(T, "all", "omitmissing");

    stat.th
    th = stat.th .* ones(size(T));
    th(self.fx + self.fy > max(self.f)) = NaN;

    clf;
    ax = gca;
    hold(ax, 'on');

    surfHandle1 = surf(ax, self.fx, self.fy, T);
    surfHandle1.FaceColor = 'interp';
    surfHandle1.EdgeColor = 'none';
    surfHandle1.FaceAlpha = 0.7;

    surfHandle2 = surf(ax, self.fx, self.fy, th);
    surfHandle2.FaceColor = 'interp';
    surfHandle2.EdgeColor = 'none';
    surfHandle2.FaceAlpha = 0.7;

    thisThreshold = th(1);
    [C, hC] = contour3(ax, self.fx, self.fy, T, [thisThreshold, thisThreshold], ...
        'LineWidth', 2, 'LineColor', 'r');
    for i = 1:numel(hC.Children)
        hC.Children(i).ZData = thisThreshold * ones(size(hC.Children(i).ZData));
    end

    xlim(ax, [min(self.f, [], 'all'), max(self.f, [], 'all')]);
    ylim(ax, [min(self.f, [], 'all'), max(self.f, [], 'all')]);
    dtick = 5;
    freqTicks = 0:dtick:max(self.f);

    view(ax, -10, 49);
    shading(ax, 'interp');
    colormap(ax, plt.viridis);
    clim(ax, [0, Tmax]);
    zlim(ax, [0, Tmax]);

    xticks(ax, freqTicks);
    yticks(ax, freqTicks);
    ax.XAxis.FontSize = FontSize;
    ax.YAxis.FontSize = FontSize;
    ax.ZAxis.FontSize = FontSize;
    box(ax, 'off');
    title(ax, '');

end