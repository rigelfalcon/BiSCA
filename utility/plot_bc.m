function plot_bc(fx, fy, bc, op, ax)
    % Plots bicoherence on a surface.
    %
    % Backward compatible signature:
    %   plot_bc(fx, fy, bc, op)
    % Tile/subplot-safe signature:
    %   plot_bc(fx, fy, bc, op, ax)
    %
    % Inputs:
    %   fx, fy: Frequency grids (2D matrices)
    %   bc: Bicoherence matrix
    %   op: Optional operator (default: @abs)
    %   ax: Optional target axes. If provided, do not call clf/gca.

    FontSize = 12;

    if nargin < 4 || isempty(op)
        op = @abs;
    end

    if nargin < 5 || isempty(ax)
        hFig = gcf;

        % Set figure position
        p = get(0, "MonitorPositions");
        if size(p, 1) > 1
            Position = p(2, :);
        else
            Position = p(1, :);
        end
        Position = [Position(1) + 200, Position(2) + 200, 550, 550];
        hFig.Position = Position;

        clf;
        ax = gca;
    else
        cla(ax);
    end

    hold(ax, 'on');
    bc = op(bc);
    f = unique(fx(:));

    surfHandle = surf(ax, fx, fy, bc);
    surfHandle.FaceColor = 'interp';
    surfHandle.EdgeColor = 'none';


    % -- 4) Configure axes --

    xlim(ax, [min(f), max(f)]);
    ylim(ax, [min(f), max(f)]);
    dtick = 5;
    freqTicks = 0:dtick:max(fx(:));

    view(ax, -10, 49);
    shading(ax, 'interp');
    colormap(ax, plt.viridis); % Use MATLAB's built-in viridis colormap
    xticks(ax, freqTicks);
    yticks(ax, freqTicks);
    ax.XAxis.FontSize = FontSize;
    ax.YAxis.FontSize = FontSize;
    ax.ZAxis.FontSize = FontSize;
    box(ax, 'off');
    title(ax, '');
    % NOTE: bc may be complex normalized bispectrum; op controls what is visualized.
    % Do not force symmetric limits here; caller should set clim() explicitly.
    % (For |beta|, limits should be [0,1]. For Re{beta}, symmetric limits make sense.)

end