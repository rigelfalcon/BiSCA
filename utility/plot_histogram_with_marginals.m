function h=plot_histogram_with_marginals(x, y,xlims,ylims,xycmp,nbin)
% Function to plot a single 2D histogram with marginal histograms
% Input: x - data vector for x-axis
%        y - data vector for y-axis

% Main plot
% figure;

if nargin<3|| isempty(xlims)
    xlims=[];
end
if nargin<4|| isempty(ylims)
    ylims=[];
end
if nargin<6|| isempty(nbin)
    nbin=30;
end



ax_main = gca;
h = histogram2(ax_main, x, y, nbin, ...
    'FaceColor', 'flat', 'EdgeColor', 'none', 'Normalization', 'pdf');
ax_main.Box = 'on';
colormap((plt.viridis))

% Optional density colorbar (avoid clutter unless explicitly enabled).
fig = ancestor(ax_main, 'figure');
if isgraphics(fig) && isappdata(fig, 'show_density_colorbar') && getappdata(fig, 'show_density_colorbar')
    cb = colorbar(ax_main);
    cb.Label.String = 'Point density (2D histogram, pdf)';
end


if ~isempty(xlims)
    xlim(xlims)
end
if ~isempty(xlims)

    ylim(ylims)
end
view(2);

% X-axis marginal histogram
ax_xMarginal = axes('Position', get_x_marginal_position(ax_main.Position));

if ~isempty(xlims)
    if isnumeric(xlims)
        histogram(ax_xMarginal, x, nbin, 'FaceColor', xycmp(1,:), 'EdgeColor', 'none','BinLimits',xlims,'Normalization','pdf');%[0.2 0.6 0.8]
    else
        histogram(ax_xMarginal, x, nbin, 'FaceColor', xycmp(1,:), 'EdgeColor', 'none','Normalization','pdf');
    end
end
xlim(xlims)
ax_xMarginal.XAxis.Visible = 'off';
ax_xMarginal.YAxis.Visible = 'off';

% Y-axis marginal histogram
ax_yMarginal = axes('Position', get_y_marginal_position(ax_main.Position));
if ~isempty(ylims)
    if isnumeric(xlims)
        histogram(ax_yMarginal, y, nbin, 'FaceColor', xycmp(2,:), 'EdgeColor', 'none', 'Orientation', 'horizontal','BinLimits',ylims,'Normalization','pdf');%[0.8 0.4 0.2]
    else
        histogram(ax_yMarginal, y, nbin, 'FaceColor', xycmp(2,:), 'EdgeColor', 'none', 'Orientation', 'horizontal','Normalization','pdf');%[0.8 0.4 0.2]
    end
end
ylim(ylims)
ax_yMarginal.XAxis.Visible = 'off';
ax_yMarginal.YAxis.Visible = 'off';
% Move marginal histograms to the bottom layer
uistack(ax_xMarginal, 'bottom');
uistack(ax_yMarginal, 'bottom');

% Set focus back to main plot
axes(ax_main);

% Improve visualization
set(gcf, 'Color', 'w');
end

function pos = get_x_marginal_position(mainPos)
    % Define the gap (adjust this value to control spacing)
    % gap = 0.004;  % gap = 0.001; Example: 2% of the figure height
    gap=0;
    % Shift the x-marginal upward by adding the gap
    pos = [mainPos(1), mainPos(2) + mainPos(4) + gap, mainPos(3), 0.04];
end


function pos = get_y_marginal_position(mainPos)
    % gap = 0.0001;  % Adjust this value to control horizontal spacing
    gap=0;
    pos = [mainPos(1) + mainPos(3) + gap, mainPos(2), 0.04, mainPos(4)];
end

% function plot_histogram_with_marginals(x, y)
%     % Function to plot a single 2D histogram with marginal histograms
%     % Input: x - data vector for x-axis
%     %        y - data vector for y-axis
%
%     % Main plot
%     histogram2(x, y, 100, 'FaceColor', 'flat');%, 'DisplayStyle', 'tile', 'ShowEmptyBins', 'off'
%     ax_main = gca;
%     ax_main.Box = 'on';
%     view(2);
%
%     % X-axis marginal histogram
%     ax_xMarginal = axes('Position', get_x_marginal_position(ax_main.Position));
%     histogram(ax_xMarginal, x, 100, 'FaceColor', 'b');
%     ax_xMarginal.XAxis.Visible = 'off';
%     ax_xMarginal.YAxis.Visible = 'off';
%
%     % Y-axis marginal histogram
%     ax_yMarginal = axes('Position', get_y_marginal_position(ax_main.Position));
%     histogram(ax_yMarginal, y, 100, 'FaceColor', 'b', 'Orientation', 'horizontal');
%     ax_yMarginal.XAxis.Visible = 'off';
%     ax_yMarginal.YAxis.Visible = 'off';
%     axes(ax_main)
%     % Improve visualization
%     set(gcf, 'Color', 'w');
% end
%
% function pos = get_x_marginal_position(mainPos)
%     % Helper function to get the position for the x-axis marginal histogram
%     pos = [mainPos(1), mainPos(2) + mainPos(4), mainPos(3), 0.05];
% end
%
% function pos = get_y_marginal_position(mainPos)
%     % Helper function to get the position for the y-axis marginal histogram
%     pos = [mainPos(1) + mainPos(3), mainPos(2), 0.05, mainPos(4)];
% end