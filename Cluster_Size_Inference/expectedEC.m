function EEC = expectedEC( S, D, detLambda, u, sigma)
% expectedEC calculates the expected Euler characteristic of an excursion
% set.
%--------------------------------------------------------------------------
% ARGUMENTS
% S     the total volume of the image.
% D     the number of dimensions.
% detLambda  the determinant of Lambda the covariance matrix of partial
%       derivatives.
% u     the threshold.
% sigma the standard derivation (assumed to be the same across the whole
%       image). Default taken to be 1.
%--------------------------------------------------------------------------
% OUTPUT
% Em    the expected value of m: the number of connected components of an 
%       excursion set.
%--------------------------------------------------------------------------
% EXAMPLES
% Dim = [250,250];
% D = length(Dim);
% S = prod(Dim);
% 
% FWHM = 10;
% detLambda = FWHM/sqrt(4*log(2)); 
% u = 2.5;
% Em = expectedm(S, D, detLambda, u)
%--------------------------------------------------------------------------
% See Alder 1981 p111-112 for the derivation of the formula.
%--------------------------------------------------------------------------
% AUTHOR: Sam Davenport.
%--------------------------------------------------------------------------
if nargin < 5
    sigma = 1;
end

if D == 1 || D == 2
    EEC = S*(2*pi)^(-(D+1)/2)*sigma^(-(2*D-1))*detLambda^(1/2)*u^(D-1)*exp(-u^2/2);
elseif D == 3
    EEC = S*(2*pi)^(-(D+1)/2)*sigma^(-(2*D-1))*detLambda^(1/2)*exp(-u^2/2)*(u^(D-1) - sigma^2);
else
    error('Only coded up for N = 1, 2, 3')
%     j = 0:((N-1)/2);
%     aj = factorial(2*j)./factorial(j)./2.^j;
%     
%     sumvec = (-ones(length(j)).^j.*aj.*sigma.^(2*j);
end

end

