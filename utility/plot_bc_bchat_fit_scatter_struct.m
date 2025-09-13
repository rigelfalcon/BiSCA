function plot_bc_bchat_fit_scatter_struct(bxa_struct, cmpath)
% Wrapper for plot_bc_bchat_fit_scatter that accepts a precomputed struct
% Usage: plot_bc_bchat_fit_scatter_struct(bixialpha)

    if nargin < 2 || isempty(cmpath)
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

    % Check required fields
    required = {'bc', 'bchat', 'fx', 'fy', 'f', 'para_alpha'};
    for i = 1:numel(required)
        if ~isfield(bxa_struct, required{i})
            error('Missing required field: %s', required{i});
        end
    end

    % Prepare data
    bc = real(bxa_struct.bc);
    bchat = real(bxa_struct.bchat);
    bc(bxa_struct.fx + bxa_struct.fy > max(bxa_struct.f)) = NaN;
    bchat(bxa_struct.fx + bxa_struct.fy > max(bxa_struct.f)) = NaN;

    clf;
    ax = gca;
    hold(ax, 'on');

    % Peak guide lines
    mu_all = bxa_struct.para_alpha.kernel.para.mu(:);
    mu_pos = sort(mu_all(mu_all > 0));
    if isempty(mu_pos)
        mu_alpha = NaN;
    else
        mu_alpha = mu_pos(1);
    end

    c_alpha = [1.00, 0.00, 0.00];
    c_h2    = [0.60, 0.00, 1.00];
    c_h3    = [1.00, 1.00, 0.00];

    for i = 1:numel(mu_all)
        mu_i = mu_all(i);
        if mu_i <= 0
            continue;
        end

        r = mu_i / mu_alpha;
        if abs(r - 1) < 0.15
            cl = c_alpha;
        elseif abs(r - 2) < 0.15
            cl = c_h2;
        elseif abs(r - 3) < 0.15
            cl = c_h3;
        else
            cl = [0.5, 0.5, 0.5];
        end

        max_f = @(x) 50 - x;
        plot([mu_i, mu_i], [0, max_f(mu_i)], 'Color', cl, 'LineWidth', 6, 'HandleVisibility', 'off');
        plot([0, max_f(mu_i)], [mu_i, mu_i], 'Color', cl, 'LineWidth', 6, 'HandleVisibility', 'off');
    end

    % Scatter for bc
    Xbc = bxa_struct.fx(:);
    Ybc = bxa_struct.fy(:);
    Zbc = bc(:);

    scatter3(ax, Xbc, Ybc, Zbc, ...
        10, Zbc, 'filled', ...
        'MarkerEdgeColor', 'none', 'DisplayName', 'Multitaper Re{beta} (signed)');

    % Surf for bchat
    surfHandle = surf(ax, bxa_struct.fx, bxa_struct.fy, bchat, 'DisplayName', 'BiSCA Re{beta} (signed)');
    surfHandle.FaceColor = 'interp';
    surfHandle.EdgeColor = 'none';
    surfHandle.FaceAlpha = 0.7;

    [C, hC] = contourf(ax, bxa_struct.fx, bxa_struct.fy, bchat, 10, ...
        'LineWidth', 1.5, 'DisplayName', 'BiSCA Re{beta} (signed)');

    xlim(ax, [min(bxa_struct.f, [], 'all'), max(bxa_struct.f, [], 'all')]);
    ylim(ax, [min(bxa_struct.f, [], 'all'), max(bxa_struct.f, [], 'all')]);

    dtick = 5;
    freqTicks = 0:dtick:max(bxa_struct.f);

    view(ax, -10, 49);
    shading(ax, 'interp');
    colormap(ax, plt.viridis);

    xticks(ax, freqTicks);
    yticks(ax, freqTicks);
    ax.XAxis.FontSize = FontSize;
    ax.YAxis.FontSize = FontSize;
    ax.ZAxis.FontSize = FontSize;
    box(ax, 'off');
    title(ax, '');

    legend

end
