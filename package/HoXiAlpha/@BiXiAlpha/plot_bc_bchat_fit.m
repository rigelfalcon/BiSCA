function plot_bc_bchat_fit(self, cmpath)

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

    bc = abs(self.bc);
    bchat = abs(self.bchat);
    bc(self.fx + self.fy > max(self.f)) = NaN;
    bchat(self.fx + self.fy > max(self.f)) = NaN;

    clf;
    ax = gca;
    hold(ax, 'on');

    N = 256;
    cmap_bc = [1, 0, 0];
    cmap_bchat = plt.viridis(N);

    z1Min = min(bc(:));
    z1Max = max(bc(:));
    z2Min = min(bchat(:));
    z2Max = max(bchat(:));
    rangeZ1 = max(z1Max - z1Min, eps);
    rangeZ2 = max(z2Max - z2Min, eps);

    normBC = (bc - z1Min) / rangeZ1;
    normBCHAT = (bchat - z2Min) / rangeZ2;

    idxBC = ones(size(bc));
    idxBCHAT = round(normBCHAT * (N - 1)) + 1;

    rgbBC = ind2rgb(idxBC, cmap_bc);
    rgbBCHAT = ind2rgb(idxBCHAT, cmap_bchat);

    surf(ax, self.fx, self.fy, bc, ...
        'CData', rgbBC, ...
        'FaceColor', 'none', ...
        'EdgeColor', 'r', ...
        'LineWidth', 1, ...
        'FaceAlpha', 1.0);

    surf(ax, self.fx, self.fy, bchat, ...
        'CData', rgbBCHAT, ...
        'FaceColor', 'texturemap', ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.7);

    xlim(ax, [min(self.f, [], 'all'), max(self.f, [], 'all')])
    ylim(ax, [min(self.f, [], 'all'), max(self.f, [], 'all')])

    dtick = 5;
    freqTicks = 0:dtick:max(self.f);

    view(ax, -7, 40);

    shading(ax, 'interp');

    xticks(ax, freqTicks);
    yticks(ax, freqTicks);
    ax.XAxis.FontSize = FontSize;
    ax.YAxis.FontSize = FontSize;
    xlabel(ax, 'Frequency (Hz)');
    ylabel(ax, 'Frequency (Hz)');

    box(ax, 'off');
    title(ax, '');

end