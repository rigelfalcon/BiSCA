function plot_potential_pca(result, low_dim_traj, num_steps)
    % PLOT_POTENTIAL_PCA - Visualizes potential landscape with proper aspect ratios
    % Inputs:
    %   result        - Struct from estimate_potential_pca
    %   low_dim_traj  - Trajectory coordinates (N x 2 matrix)

    % Extract data from result struct
    X1 = result.X1;
    X2 = result.X2;
    U = result.U;

    % Create figure
    fig = gcf;
    set(fig, 'DefaultAxesFontName', 'Arial', 'DefaultAxesFontSize', 12);
    ax = gca;
    % Set aspect ratio before plotting
    data_range = [range(X1(:)), range(X2(:)), range(U(:))];

    % Surface plot
    h_surf = surf(X1, X2, U, 'EdgeColor', 'none');
    shading interp
    hold on;
    % Axis styling
    view(-45, 30);
    xlim([min(X1(:)), max(X1(:))])
    ylim([min(X2(:)), max(X2(:))])

    grid off;
    box off;
    FontSize = 12;
    % Labels and colorbar
    xlabel('PC1', 'FontSize', FontSize) %, 'FontWeight', 'bold'
    ylabel('PC2', 'FontSize', FontSize)
    zlabel('Potential Landscape', 'FontSize', FontSize)

    % Add trajectory with arrow3
    if nargin > 1 && ~isempty(low_dim_traj)
        if nargin > 2 && ~isempty(num_steps)
        else
            num_steps = 2 * result.num_steps;
        end
        low_dim_traj = low_dim_traj(1:num_steps, :);
        % Upsample trajectory
        Fs_original = 200;
        Fs_upsampled = 4000;
        low_dim_traj = reshape(resampleGridded(reshape(low_dim_traj, size(low_dim_traj, 1), []), Fs_original, Fs_upsampled), [], size(low_dim_traj, 2), size(low_dim_traj, 3));


        % Get potential values along trajectory
        F = griddedInterpolant(X1', X2', U');
        U_traj = F(low_dim_traj(:, 1), low_dim_traj(:, 2));
        U_traj = U_traj + 0.1 * data_range(3);
        % Convert to 3D coordinates
        traj_3d = [low_dim_traj, U_traj];


        % Plot trajectory segments with arrow3
        arrow_scale = 0.15 * min(data_range);
        plot3(traj_3d(:, 1), traj_3d(:, 2), traj_3d(:, 3), ...
            'Color', [1 0 0 0.8], 'LineWidth', 2); % Red with 50% transparency

        interval_arrow = 120;
        interval_inner = 10;

        pos_start = traj_3d(20 + 1:interval_arrow:end - interval_arrow + 1, :);
        pos_end = traj_3d(20 + 1 + interval_inner:interval_arrow:end - interval_arrow + 1 + interval_inner, :);
        valid = sum(abs(pos_end - pos_start).^2, 2) > 0;
        pos_start = pos_start(valid, :);
        pos_end = pos_end(valid, :);

        pos_start = pos_start(1:min(5, size(pos_start, 1)), :);
        pos_end = pos_end(1:min(5, size(pos_end, 1)), :);

        arrow3(gather(pos_start), gather(pos_end), ...
            'r1r');
        % Plot start and end markers
        plot3(traj_3d(1, 1), traj_3d(1, 2), traj_3d(1, 3), ...
            'o', 'MarkerSize', 10, 'LineWidth', 2, ...
            'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
        plot3(traj_3d(end, 1), traj_3d(end, 2), traj_3d(end, 3), ...
            's', 'MarkerSize', 10, 'LineWidth', 2, ...
            'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');

    end

    colormap(redbluecmap)

    c = colorbar;
    c.Label.String = 'Potential Landscape';
    c.Label.FontSize = 16;


    hold off;
end