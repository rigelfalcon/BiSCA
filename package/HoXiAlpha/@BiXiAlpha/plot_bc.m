function plot_bc(self, comp, op, cmpath)

    if nargin < 4 || isempty(cmpath)
        cmpath = [];
    end
    hFig = gcf;

    FontSize = 6;
    LineWidth = 1;
    p = get(0, "MonitorPositions");
    if size(p, 1) > 1
        Position = p(2, :);
    else
        Position = p(1, :);
    end
    Position = [Position(1) + 200, Position(2) + 200, 160, 160];

    hFig.Position = Position;

    bc = op(self.(comp));
    bc(self.fx + self.fy > max(self.f)) = NaN;
    surfbc(self.fx, self.fy, bc);

    if sum(bc > 0, "all") && sum(bc < 0, "all")
        clim([-0.6, 0.6]);
        if isempty(cmpath)
            colormap(flip(plt.viridis))
        end
    else
        clim([0.05, 0.6]);
        cbh = colorbar('h');
        get(get(cbh, 'Children'))
        if isempty(cmpath)
            colormap(plt.viridis)
        end
    end

    xlim([min(self.f, [], 'all'), max(self.f, [], 'all')])
    ylim([min(self.f, [], 'all'), max(self.f, [], 'all')])

    dtick = 5;
    range = [0:dtick:max(self.f)];
    hold on;
    plot_bsgrid(range)
    title([])
    colorbar off
    shading interp
    xticks(range)
    yticks(range)
    hA = gca;
    hA.XAxis.FontSize = FontSize;
    hA.YAxis.FontSize = FontSize;
    xlabel('Frequency (Hz)')
    ylabel('Frequency (Hz)')

end