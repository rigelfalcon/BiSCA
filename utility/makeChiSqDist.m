function pd = makeChiSqDist(nu)
%MAKECHISQDIST Create a chi-square distribution using Gamma parameters.
%   PD = MAKECHISQDIST(NU) returns a GammaDistribution object configured
%   to represent a chi-square distribution with NU degrees of freedom.
%
%   Input:
%       nu - Degrees of freedom (positive scalar)
%
%   Output:
%       pd - GammaDistribution object parameterized as chi-square(nu)
%
%   Example:
%       pd = makeChiSqDist(5);
%       x = 0:0.1:15;
%       plot(x, pd.pdf(x));
%       title('Chi-Square PDF (ν = 5) via Gamma');

% Validate input
if ~isscalar(nu) || nu <= 0
    error('makeChiSqDist:InvalidInput', ...
        'nu must be a positive scalar.');
end

% Parameterize Gamma to represent chi-square(nu)
a = nu / 2;   % Gamma shape parameter
b = 2;        % Gamma scale parameter

% Create Gamma distribution object
pd = makedist('Gamma', 'a', a, 'b', b);

end
