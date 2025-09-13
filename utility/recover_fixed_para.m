function [x, lb, ub] = recover_fixed_para(xsub, lbsub, ubsub, parainfo)

    x = parainfo.x;
    x(~parainfo.xIndices.fixed) = xsub;

    if ~isempty(lbsub)
        lb = parainfo.lb;
        lb(~parainfo.xIndices.fixed) = lbsub;
    else
        lb = [];
    end
    if ~isempty(ubsub)
        ub = parainfo.ub;
        ub(~parainfo.xIndices.fixed) = ubsub;
    else
        ub = [];
    end

end