function [Y] = gen_ar_ts(A, Y, RES)
    %% VAR (vector autoregressive model)
    % model:
    % Y_t=c+A_1*Y_t-1+A_2*Y_t-2+...A_p*Y_t-p+t
    % A: n*n*p
    % Y: t*n*m
    % RES: t*n*m

    [n, ~, p] = size(A);
    [t, ~, m] = size(Y);
    A = reshape(flip(A, 3), [n, n * p]); %1:p is t-1:t-p
    Y = reshape(permute(Y, [2, 1, 3]), n * t, m); %t,n,m->n*t,m
    RES = permute(RES, [2, 3, 1]); %t,n,m->n,m,t

    for i = 1:t
        if i - p <= 0
            continue
        end
        % n*(i-1)-n*p -> n*(i-1)
        Y(n * (i - 1) + 1:n * i, :) = A * Y(n * (i - 1) - n * p + 1:n * (i - 1), :) + RES(:, :, i);
    end
    Y = permute(reshape(Y, n, t, m), [2, 1, 3]); %n*t,m->t,n,m

end