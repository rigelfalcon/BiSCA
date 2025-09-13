function Jinv = pinv_ridge(J, lambda)
    N = size(J, 1);
    lambda = lambda * sum(J .* J, 'all') / N;
    Jinv = (J' * J + lambda * eye(size(J, 2))) \ J';
end