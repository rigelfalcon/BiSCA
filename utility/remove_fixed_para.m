function [x, lb, ub, parainfo] = remove_fixed_para(x, lb, ub)
    nvar = length(x);
    if nvar ~= length(lb) || nvar ~= length(ub)
        error('x, lb, and ub must have the same length');
    end
    parainfo.xIndices = classifyBoundsOnVars(lb, ub, nvar, true);
    parainfo.x = x;
    parainfo.lb = lb;
    parainfo.ub = ub;

    x(parainfo.xIndices.fixed) = lb(parainfo.xIndices.fixed);

    x = x(~parainfo.xIndices.fixed);
    lb = lb(~parainfo.xIndices.fixed);
    ub = ub(~parainfo.xIndices.fixed);

end