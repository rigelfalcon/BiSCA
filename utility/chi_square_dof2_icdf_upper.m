function x = chi_square_dof2_icdf_upper(p)
    % chi_square_icdf computes the ICDF for Chi-Square distribution with dof = 2
    % Input:
    %   p : Cumulative probability (between 0 and 1)
    % Output:
    %   x : Value corresponding to the given cumulative probability p
    
    % Check if input p is between 0 and 1
    if p < 0 || p > 1 %<= >=
        error('Input p must be between 0 and 1');
    end
    
    % Calculate the ICDF using the formula for Chi-Square with dof = 2
    % x = -2 * log(1 - p);
    % syms x
    % f=taylor(log(x), x, 'ExpansionPoint', 1, 'Order', 1000);
    % f=matlabFunction(f);
    % x=-2*f(1-p)
    x = -2 * log(p);
end
