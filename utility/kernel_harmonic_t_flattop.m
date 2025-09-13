function [y, ycomp] = kernel_harmonic_t_flattop(x, h, mu, sigma, nu, d, b)
    p = 2;

    if size(sigma, 2) == 3
        x = [x, sum(x, 2)]; %x,y,x+y
        mu = [mu, sum(mu, 2)];
        p = 1;
    end

    sigma = permute(sigma, [3, 1, 2]); % Nmu*Nd*1->1*Nmu*Nd
    mu = permute(mu, [3, 1, 2]); % Nmu*Nd*1->1*Nmu*Nd
    nu = permute(nu, [3, 1, 2]); % Nmu*Nd*1->1*Nmu*Nd
    d = permute(d, [3, 1, 2]); % Nmu*Nd*1->1*Nmu*Nd
    b = permute(b, [3, 1, 2]); % Nmu*Nd*1->1*Nmu*Nd

    x = permute(x, [1, 3, 2]); % Nx*Nd->Nx*1*Nd

    if nargout > 1
        ycomp = h .* (prod((1 + abs((x - mu) ./ sigma) .^ d) .^ (-nu), 3) .^ p);
        y = sum(ycomp, 2);

    else
        y = sum(h .* prod((1 + abs((x - mu) ./ sigma) .^ d) .^ (-nu), 3) .^ p, 2);

    end

end

function sigma = scale_sigma(sigma, h, nu, d, p)
    % decouple sigma and nu
    sigma = sigma ./ ((2) .^ (1 ./ (nu)) - 1) .^ (1 ./ d);

end

function y = tkernel(x, mu, sigma, nu, d)
    y = (1 + abs((x - mu) ./ sigma) .^ d) .^ (-1 ./ nu);
end