function inactiveVars = find_inactive_variables(fun, x0, epsilon, tolerance)
% Determines inactive variables for a function at a point x0.
% fun - Function handle, should accept a vector x of the same size as x0.
% x0 - Vector of parameters at which the function is evaluated.
% epsilon - Small perturbation applied to each parameter to test sensitivity.
% tolerance - Tolerance within which a change is considered negligible.
    if nargin < 3
        epsilon = 1e-6;  % Default perturbation
    end
    if nargin < 4
        tolerance = 1e-6; % Default tolerance
    end
    n = numel(x0);
    inactiveVars = true(1, n);  % Assume all variables are inactive initially
    f0 = fun(x0);               % Evaluate function at the base point

    for i = 1:n
        xTest = x0;
        xTest(i) = xTest(i) + epsilon; % Perturb each variable slightly
        fTest = fun(xTest);            % Evaluate function at the perturbed point
        % Check if the change in function output is negligible
        if any(abs(fTest - f0) > tolerance)
            inactiveVars(i) = false;  % Variable is active
        end
    end
end
