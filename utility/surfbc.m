function surfbc(fx, fy, bc)

p = surf(fx, fy, bc);
view(0, 90)
xlim([min(fx, [], 'all'), max(fx, [], 'all')])
ylim([min(fy, [], 'all'), max(fy, [], 'all')])
set(p, 'EdgeColor', 'none')
% set(p, 'EdgeColor', 'none', 'FaceAlpha', 0.5)

% shading interp
grid off
ax = gca;
% axis equal
% colormap(turbo)% colormap(flip(hot))%hot
% colormap(ax, flip(hot))
colormap(ax, (hot))

% colorbar
% title('bicoherence')
box off
% clim([0.05, 0.6])
clim([0.05, 0.8])

end
