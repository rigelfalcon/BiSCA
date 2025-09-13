function plot_bc_bchat_fit_scatter(self, cmpath)

    % -- 1) Handle inputs & basic figure setup --
    if nargin < 4 || isempty(cmpath)
        cmpath = [];
    end

    hFig = gcf;
    FontSize = 16;
    LineWidth = 1;

    % Optionally set the figure position
    p = get(0, "MonitorPositions");
    if size(p, 1) > 1
        Position = p(2, :);
    else
        Position = p(1, :);
    end
    Position = [Position(1) + 200, Position(2) + 200, 550, 550];
    hFig.Position = Position;

    % -- 2) Prepare data: absolute value, then NaN-out region above max(f) --
    bc = real(self.bc);
    bchat = real(self.bchat);
    bc(self.fx + self.fy > max(self.f)) = NaN;
    bchat(self.fx + self.fy > max(self.f)) = NaN;

    % -- 3) Create a single axes and "hold on" so both plots appear together --
    clf;
    ax = gca;
    hold(ax, 'on');

    %%
% Colorblind-safe palette for peak guide lines 
% Use colors that avoid conflict with viridis colormap.
mu_all = self.para_alpha.kernel.para.mu(:);
mu_pos = sort(mu_all(mu_all > 0));
if isempty(mu_pos)
    mu_alpha = NaN;
else
    mu_alpha = mu_pos(1);
end

% Peak colors - avoid inferno/viridis overlap
c_alpha = [1.00, 0.00, 0.00];  % alpha peak -> Pure Red
c_h2    = [0.60, 0.00, 1.00];  % 2/1 harmonic -> Bright Violet
c_h3    = [1.00, 1.00, 0.00];  % 3/1 harmonic -> Bright Yellow

for i = 1:numel(mu_all)
    mu_i = mu_all(i);
    if mu_i <= 0
        continue;
    end

    % Classify by harmonic ratio relative to the smallest positive mu.
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

    % -- 4) Scatter3 for bc --
    Xbc = self.fx(:);
    Ybc = self.fy(:);
    Zbc = bc(:);

    scatter3(ax, Xbc, Ybc, Zbc, ...
        10, ... % Marker size
        Zbc, ... % Color each point by bc value
        'filled', ... % Solid markers
        'MarkerEdgeColor', 'none', 'DisplayName', 'Multitaper Re{beta} (signed)'); %'none'

    % -- 5) Surf for bchat --
    surfHandle = surf(ax, self.fx, self.fy, bchat, 'DisplayName', 'BiSCA Re{beta} (signed)');

    surfHandle.FaceColor = 'interp';
    surfHandle.EdgeColor = 'none';
    surfHandle.FaceAlpha = 0.7; % Slight transparency for overlap

    [C, hC] = contourf(ax, self.fx, self.fy, bchat, 10, ...
        'LineWidth', 1.5, 'DisplayName', 'BiSCA Re{beta} (signed)');

    % -- 6) Axes limits, grid, etc. --
    % Shared frequency limits
    xlim(ax, [min(self.f, [], 'all'), max(self.f, [], 'all')]);
    ylim(ax, [min(self.f, [], 'all'), max(self.f, [], 'all')]);

    % Example: custom ticks if desired
    dtick = 5;
    freqTicks = 0:dtick:max(self.f);

    % 3D angled view (change to (0,90) for top-down)
    view(ax, -10, 49);

    % Interpolated shading
    shading(ax, 'interp');

    % Single colormap & colorbar for both objects
    colormap(ax, plt.viridis);

    % Configure ticks and fonts
    xticks(ax, freqTicks);
    yticks(ax, freqTicks);
    ax.XAxis.FontSize = FontSize;
    ax.YAxis.FontSize = FontSize;
    ax.ZAxis.FontSize = FontSize;
    box(ax, 'off');
    title(ax, '');

    legend


end