function x_estimate = gdescent( initial_estimate, fderiv, gamma, tol, f, max_iters)
% NEWFUN serves as a function template.
%--------------------------------------------------------------------------
% ARGUMENTS
% 
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
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
    x_update = x_estimate - gamma*fderiv(x_estimate);
    delta = norm(fderiv(x_estimate));
    x_estimate = x_update;
    if delta < tol
        break
    end
end

end
