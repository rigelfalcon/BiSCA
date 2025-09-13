function [shat, scomp] = pred_s_alpha(f, ks, kh, para_alpha)
    % scomp: different with comp, scomp is with the scale factor, but comp is
    %        unit height
    if nargout > 1
        [shat, scomp] = para_alpha.kernel.eval(f);
    else
        [shat] = para_alpha.kernel.eval(f);
    end
end