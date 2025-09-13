function plot_para_dist( pd, tval, hval)

if isfield(pd,'DistributionName')&&strcmpi(pd.DistributionName, 'Uniform')
    grid off
    box off
    return
end

hold on;

if hval
    h = pd.plot;
    h(2).FaceColor = 'red';
    xline(tval, 'LineWidth', 2, 'Color', 'blue')

else
    h = pd.plot;
    h(2).FaceColor = 'green';
    xline(tval, 'LineWidth', 2, 'Color', 'blue')
end

hold off;
end