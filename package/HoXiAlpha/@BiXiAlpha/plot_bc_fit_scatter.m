function plot_bc_fit_scatter(self, cmpath)

    if nargin < 4 || isempty(cmpath)
        cmpath = [];
    end

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

    % NOTE: self.bc is complex normalized bispectrum (complex bicoherence).
    % We plot its real part, so negative values are expected.
    bc = real(self.bc);
    bc(self.fx + self.fy > max(self.f)) = NaN;

    clf;
    ax = gca;
    hold(ax, 'on');

    surfHandle = surf(ax, self.fx, self.fy, bc);
    surfHandle.FaceColor = 'interp';
    surfHandle.EdgeColor = 'none';
    surfHandle.FaceAlpha = 0.7;

    [C, hC] = contourf(ax, self.fx, self.fy, bc, 10, ...
        'LineWidth', 1.5);

    xlim(ax, [min(self.f, [], 'all'), max(self.f, [], 'all')]);
    ylim(ax, [min(self.f, [], 'all'), max(self.f, [], 'all')]);

    dtick = 5;
    freqTicks = 0:dtick:max(self.f);

    view(ax, -10, 49);

    shading(ax, 'interp');

    colormap(ax, plt.viridis);

    % NOTE (Fig.1D): Legend for the colored center lines (alpha peak / 2/1 / 3/1)
    % is handled in the manuscript assembly (see doc/reply/archive/R3_Response_with_Appendix.md).

    xticks(ax, freqTicks);
    yticks(ax, freqTicks);
    ax.XAxis.FontSize = FontSize;
    ax.YAxis.FontSize = FontSize;
    ax.ZAxis.FontSize = FontSize;
    box(ax, 'off');
    title(ax, '');

end