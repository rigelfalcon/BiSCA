function x0 = strict_bound_x0(x0, lb, ub)
    x0 = x0(:);
    lb = lb(:);
    ub = ub(:);

    if length(x0) ~= length(lb) || length(x0) ~= length(ub)
        error('x0, lb, and ub must have the same length.')
    end

    if any(lb > ub)
        error('Lower bounds must be less than or equal to upper bounds.')
    end

    nvar = length(x0);
    xIndices = classifyBoundsOnVars(lb(:), ub(:), nvar, true);
    violatedFixedBnds_idx = x0(xIndices.fixed) ~= lb(xIndices.fixed);
    violatedLowerBnds_idx = x0(xIndices.finiteLb) <= lb(xIndices.finiteLb);
    violatedUpperBnds_idx = x0(xIndices.finiteUb) >= ub(xIndices.finiteUb);
    xIndices.idFixed = find(xIndices.fixed);
    xIndices.idFiniteLb = find(xIndices.finiteLb);
    xIndices.idFiniteUb = find(xIndices.finiteUb);
    violatedFixedBnds_idx = xIndices.idFixed(violatedFixedBnds_idx);
    violatedLowerBnds_idx = xIndices.idFiniteLb(violatedLowerBnds_idx);
    violatedUpperBnds_idx = xIndices.idFiniteUb(violatedUpperBnds_idx);

    if any(violatedFixedBnds_idx)
        x0(violatedFixedBnds_idx) = lb(violatedFixedBnds_idx);
    end
    if any(violatedLowerBnds_idx)
        x0(violatedLowerBnds_idx) = lb(violatedLowerBnds_idx) + (1e-4) * (ub(violatedLowerBnds_idx) - lb(violatedLowerBnds_idx));
    end
    if any(violatedUpperBnds_idx)
        x0(violatedUpperBnds_idx) = ub(violatedUpperBnds_idx) - (1e-4) * (ub(violatedUpperBnds_idx) - lb(violatedUpperBnds_idx));
    end

end