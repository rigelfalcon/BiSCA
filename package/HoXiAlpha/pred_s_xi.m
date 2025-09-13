function [shat, scomp] = pred_s_xi(f, para_xi)
    if nargout > 1
        [shat, scomp] = para_xi.kernel.eval(f);
    else
        [shat] = para_xi.kernel.eval(f);
        scomp = [];
    end
end