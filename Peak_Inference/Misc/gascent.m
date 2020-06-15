function x_estimate = gascent( initial_estimate, fderiv, gamma, tol, f, max_iters)
% GASCENT( initial_estimate, fderiv, gamma, tol, f, max_iters) performs
% gradient ascent from an initial starting point and ends at a local
% maximum.
%--------------------------------------------------------------------------
% ARGUMENTS
% initial_estimate  the starting point for the search
% fderiv   	a function handle giving the derivative of the function
%          	you wish to maximize
% gamma 	the (for now fixed) parameter used to perform gradient ascent
% tol       the size of the tolerance. The algorithm stops once
%           min(norm(x_update - x_estimate,2), norm(fderiv(x_estimate))) < tol
% f         the original function (optional argument to check progress)
% max_iters     the maximum number of iterations performed until the algorithm
%               will terminate. (Default is 1000.)
%--------------------------------------------------------------------------
% OUTPUT
% x_estimate    the estimate of the location of the local maximum
%--------------------------------------------------------------------------
% EXAMPLES
% gascent( 0.1, @(x)-2*x )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
if nargin < 4
    gamma = 0.2;
end
if nargin < 5
    tol = 10^(-4);
end
if nargin < 7
    max_iters = 1000;
end

x_estimate = initial_estimate;

for I =1:max_iters
    x_update = x_estimate + gamma*fderiv(x_estimate); %NOTE we're doing ascent so want a plus here!!
    %     fderiv(x_estimate)
%     f(x_update)
    %     delta = norm(x_update - x_estimate,2) + norm(fderiv(x_estimate));
    if nargin < 5
        delta = min(norm(x_update - x_estimate,2), norm(fderiv(x_estimate)));
    else
        delta = min(min(norm(x_update - x_estimate,2), norm(fderiv(x_estimate))), abs(f(x_update) - f(x_estimate)));
    end
%     delta = min(norm(x_update - x_estimate,2), norm(fderiv(x_estimate)));
%     delta = min(norm(fderiv(x_estimate)), norm(f(x_estimate) - f(x_update)));
%     delta = max(norm(x_update - x_estimate,2), norm(f(x_estimate) - f(x_update)));
    x_estimate = x_update;
    if delta < tol
        break
    end
end

end
