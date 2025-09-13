function surfbs(fx, fy, bs, issymcolor)
    if nargin < 4 || isempty(issymcolor)
        issymcolor = true;
    end
    p = surf(fx, fy, bs);
    view(0, 90)
    xlim([min(fx, [], 'all'), max(fx, [], 'all')])
    ylim([min(fy, [], 'all'), max(fy, [], 'all')])
    set(p, 'edgecolor', 'none')
    grid off
    ax = gca;
    colormap(flip(hot))

    colorbar
    if issymcolor
        maxbs = max(abs(bs), [], "all");
    end
    title('bispectrum')
    box off
end