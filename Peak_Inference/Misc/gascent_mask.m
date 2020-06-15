function x_estimate = gascent_mask( initial_estimate, fderiv, mask, gamma, tol, f, max_iters)
% GASCENT_MASK( initial_estimate, fderiv, gamma, tol, f, max_iters) performs
% gradient ascent within a mask from an initial starting point and ends at 
% a local maximum (within the mask). Note it assumes that the lattice
% points of the msak are evenly space integers.
%--------------------------------------------------------------------------
% ARGUMENTS
% initial_estimate  the starting point for the search
% fderiv   	a function handle giving the derivative of the function
%          	you wish to maximize
% mask      a 0-1 array with 1 indicating a voxel lieing within the image
%           and 0 indicating a voxel outside the mask. 
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
% Y = [1,1,1;1,2,1;1,1,1]
% mask = [0,0,1;0,0,1;0,0,1];
% FWHM = 3;
% Kprime = @(x) GkerMVderiv(x,FWHM);
% field_deriv = @(tval) applyconvfield(tval, Y, Kprime, 0 );
% initial_estimate = [3,3]';
% gascent(initial_estimate, field_deriv) %Gets to 2 as expected!
% gascent_mask(initial_estimate, field_deriv, mask) %nails it perfectly!
% f = @(tval) applyconvfield(tval, Y, FWHM, 0 );
% gascent_mask(initial_estimate, field_deriv, mask, 10)
% % note that if you look at how the x_updates changes, it's trying to send
% % the point further towards the maximum and then projects it back!
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
if nargin < 4
    gamma = 0.1;
end
if nargin < 5
    tol = 10^(-4);
end
if nargin < 7
    max_iters = 1000;
end

x_estimate = initial_estimate;

for I = 1:max_iters
%     x_estimate;
    x_update = x_estimate + gamma*fderiv(x_estimate);
    if norm(x_update - x_estimate) > 0.5
       x_update = (x_update - x_estimate)/(2*norm(x_update - x_estimate)) + x_estimate;
    end
%     pause
%     x_update
%     x_update
    if ~inmask(x_update, mask)
        x_update = maskprojection(x_update, mask);
    end
%     f(x_update)
%     if x_update(1) > 10.5 
%         a = 4
%     end
    if nargin < 6
        delta = min(min(norm(x_update - x_estimate,2), norm(fderiv(x_estimate))), abs(f(x_update) - f(x_estimate)) );
    else
        delta = min(min(norm(x_update - x_estimate,2), norm(fderiv(x_estimate))), abs(f(x_update) - f(x_estimate)) );
    end
    x_estimate = x_update;
%     if nargin > 6
%         f(x_estimate)
%     end
    if delta < tol
        break
    end
end

end
