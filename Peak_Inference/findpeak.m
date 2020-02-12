function x_estimate = findpeak( initial_guess, fprime, fprime2, mask, max_or_min, tol, gamma, g_ascent_tol, f )
% FINDPEAK( initial_guess, f, fprime, fprime2, max_or_min, tol, max_iters )
% finds the peak (either a maximum or minimum) of a random field by
% searching from an initial estimate.
%
% Basically we gradient ascend our way out of minima and saddle points and
% newton raphson our way into maxima. (As the gradient descent algorithm
% gets slower as you approach a maximum and NR can get stuck at local
% minima or saddle points!) Plus near saddle points/local minima you should 
% be taking the line of steepest ascent up. Near boundaries should just use
% gradient ascent as the derivative is not zeor so you won't converge
% slowly!! May need to investigate peturbing things if you still get stuck
% in saddle points.  NEED TO Put In A CLAUSE FOR WHEN THE FPRIME2 is not
% specfied forcingit to use gradient ascent in this case!.
%--------------------------------------------------------------------------
% ARGUMENTS
% initial_guess     an initial guess of the peak location
% f                 a function handle representing the function
% fprime            a function handle representing the function's derivative
%                   note that if fprime is NaN it is estimated using f.
% fprime2           a function handle representing the function's 2nd derivative
%                   note that if fprime2 is NaN it is estimated using fprime
% max_or_min        0/1 specifies whether to search for a local maximum or
%                   minimum. Default is 1: local maximum. (NEEDS MORE
%                   IMPLEMENATION FOR MINiMA!)
% tol               a tolerance representing the difference of the function
%                   from 0
% max_iters         the maximum number of iterations that the algorithm can take
%--------------------------------------------------------------------------
% OUTPUT
% x_estimate        an estimate of the location of the zero of the function
%--------------------------------------------------------------------------
% EXAMPLES
% R -> R
% findpeak(0.1, @(x) x^4, @(x)4*x^3, @(x) 12*x^2, 1)
%
% Y = [1,2,1;1,1,1;1,1,1] %I.e so the peak will be outside the mask!
% mask = [0,1,1;0,1,1;0,1,1];
% FWHM = 3;
% Kprime = @(x) GkerMVderiv(x,FWHM);
% Kprime2 = @(x) GkerMVderiv2(x,FWHM);
% field = @(tval) applyconvfield(tval, Y, FWHM, 0 );
% field_deriv = @(tval) applyconvfield(tval, Y, Kprime, 0 );
% D = 2;
% field_deriv2 = @(tval) reshape(applyconvfield(tval, Y, Kprime2), [D,D]);
% initial_estimate = [3,3]';
% peakloc = findpeak( initial_estimate, field_deriv, field_deriv2, mask )
% field([1,2]')
% field(peakloc)
% field(peakloc +0.01)
% field(peakloc -0.01)
%
% Y = [1,1,1;10,1,1;1,1,1] %I.e so the peak will be outside the mask!
% field = @(tval) applyconvfield(tval, Y, FWHM, 0 );
% field_deriv = @(tval) applyconvfield(tval, Y, Kprime, 0 );
% D = 2;
% field_deriv2 = @(tval) reshape(applyconvfield(tval, Y, Kprime2), [D,D]);
% peakloc = findpeak( [3,3]', field_deriv, field_deriv2, mask )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
if nargin < 3
    fprime2 = NaN;
end
if nargin < 4
    mask = [];
end
if nargin < 5
    max_or_min = 1;
end
if nargin < 6
    tol = 10^(-5);
end
if nargin < 7
    gamma = 0.05;
end
if nargin < 8
    g_ascent_tol = 0.01;
end
if max_or_min == 0 %I.e. allow 0 to be used to denote minima
    max_or_min = -1;
end

if isempty(mask)
    if max_or_min
        ga_estimate = gascent( initial_guess, fprime, gamma, 0.1); %high valued tol
    else
        ga_estimate = gdescent( initial_guess, fprime, gamma, 0.1); %high valued tol
    end
    x_estimate = NewtonRaphson( fprime, ga_estimate, fprime2, tol );
else
    %For now we're assuming the mask is ones everywhere on a rectangle so
    %the boundary is just the outer rectangle! Need to generalize this of 
    %course.
    D = length(size(mask));
    boundary = zeros(size(mask)); 
    if D == 1
        boundary(1) = 1; boundary(end) = 1;
    elseif D == 2
        boundary(1, :) = 1; boundary(end, :) = 1; boundary(:,1) = 1; boundary(:, end) = 1;
    end
    if inmask(initial_guess, boundary)
        x_estimate = gascent_mask(initial_guess, fprime, mask, gamma, 0.00001, f);
        % maybe should tkae a large gamma but enforce the algorithm to go
        % no more than 0.5 voxels each time!
    else
        ga_estimate = gascent( initial_guess, fprime, gamma, g_ascent_tol, f); %high valued tol
%         ga_estimate = gascent( initial_guess, fprime, 0.01, 0.00001, f); %high valued tol
        x_estimate = NewtonRaphson( fprime, ga_estimate, fprime2, tol );
    end
end

% NEED TO IMPLEmENT The below to CATCH IF THE END POINT NOT ON THE MASK
% however don't need to run it if you tried the mask above. Trouble is
% otherwise you don't have a mask so the below stuff doesn't quite work!
% Fix it.
% if ~inmask(x_estimate, mask) %I.e if the above failed to find a point that actually lies in the mask take action!
%     x_estimate = gascent_mask(initial_guess, fprime, mask, 0.05, tol);
% end

end

% while difference > tol
%     x_estimate
%     f(x_estimate)
%     gradmate = fprime2(x_estimate);
%     [~,D] = eig(gradmate)
%     any(max_or_min*diag(D) >= 0)
%     pause
% %     D
%     if any(max_or_min*diag(D) >= 0)
%         % Can do more than 1 round of gradient ascent here!
%         gamma = 0.1;
%         x_estimate = x_estimate + gamma*fprime(x_estimate);
%     else
%         h = -pinv(gradmate)*fprime(x_estimate);
%         x_estimate = x_estimate + h;
%     end
%     
%     difference = sum(abs(fprime(x_estimate)));
%     % Note that this stopping criterion means that if you start at a minima
%     % you'll stay there even if you want to get to a maxima. But hopefully
%     % you can be clever enough not to initialize the algorithm there!!
%     
%     iters = iters + 1;
%     if iters > max_iters
%         x_estimate = NaN*ones(1, D);
%         return
%     end
% end