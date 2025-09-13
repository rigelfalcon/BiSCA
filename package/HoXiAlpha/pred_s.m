function [shat, scomp] = pred_s(f, ks, kh, para_xi, para_alpha)
    % scomp: different with comp, scomp is with the scale factor, but comp is
    %        unit height

    if nargout > 1
        [s_xi, scomp_xi] = pred_s_xi(f, para_xi);

        [s_alpha, scomp_alpha] = pred_s_alpha(f, ks, kh, para_alpha); %,comp gaussian_kernel(f,mu,sigma).*para_alpha.kernel.alpha.';
        shat = s_xi(:) + s_alpha;
        scomp = [scomp_xi, scomp_alpha];
    else
        [s_xi] = pred_s_xi(f, para_xi);
        [s_alpha] = pred_s_alpha(f, ks, kh, para_alpha); %,comp gaussian_kernel(f,mu,sigma).*para_alpha.kernel.alpha.';
        shat = s_xi(:) + s_alpha;
        scomp = [];
    end
end