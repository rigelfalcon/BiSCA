function c = cumulants( x, n )
    % this function returns the cumulants of the input signal x up to order n
    % c = cumulants( x, n )
    % x: input signal
    % n: order of cumulants
    % c: cumulants of x up to order n

    % check input
    if nargin < 2
        error('Not enough input arguments.');
    end
    if n < 1
        error('Order of cumulants must be greater than 0.');
    end
    
    % compute cumulants
    c = zeros(n,size(x,2));
    for i = 1:n
        c(i,:) = mean( x.^i );
    end
end
