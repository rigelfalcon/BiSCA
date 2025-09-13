function [ci, se] = nlpar_ci_se_bc(bc, resid, jac_r2x, jac_bc2x)
    lambda = 1e-3;
    alpha = 0.05;

    Nr = size(jac_r2x, 1);
    Nbc = size(jac_bc2x, 1);
    v = Nr - Nbc;
    rmse = norm(resid) .^ 2 / v;
    jac_r2x_inv = pinv_ridge(jac_r2x, lambda);
    se = sqrt(rmse * diag(jac_bc2x * (jac_r2x_inv * jac_r2x_inv') * jac_bc2x'));

    % Calculate confidence interval
    delta = se * tinv(1 - alpha / 2, v);
    ci = [(bc(:) - delta) (bc(:) + delta)];
end