function x_estimate = NewtonRaphson( f, initial_guess, fprime, tol, max_iters )
% NewtonRaphson finds the zeros of a uni or multivariate function given an
% initial guess and the form of its derivative.
%--------------------------------------------------------------------------
% ARGUMENTS
% f                 a function handle representing the function
% initial_guess     an initial guess of the zero location
% fprime            a function handle representing the function's derivative
%                   note that if fprime is NaN it is estimated using f.
% tol               a tolerance representing the difference of the function
%                   from 0
% max_iters         the maximum number of iterations that the algorithm can take
%--------------------------------------------------------------------------
% OUTPUT
% x_estimate        an estimate of the location of the zero of the function
%--------------------------------------------------------------------------
% EXAMPLES
% R -> R
% NewtonRaphson(@(x) x^2 - 612, @(x)2*x, 1)
% NewtonRaphson(@examplefunction, @examplefunctionprime, 1)
% NewtonRaphson(@(theta) loglikTNderiv(2,theta,1), @(theta) loglikTNderivprime(2,theta,1), 2)
%
% %R -> RN No point doing this one as may as well zero the first coordinate
% %using univariate Newton Raphson.
%
% %R2 -> R
% NewtonRaphson(@(x) x(1) + x(2), [1,2], @(x) [1,1])
%
% Note this version is quite robust even when the derivative is zero at the
% zero!
% Taylor's theorem e.g.:
% NewtonRaphson(@(x) x(1)^2 + x(2)^2, [0,1], @(x)[2*x(1),2*x(2)], 10^(-6))
%
% Multidimensinal output and input:
% NewtonRaphson(@(x) [x(1)^2 - 1; x(1)^2 + x(2)^2], [10,5], @(x)[2*x(1),0; 2*x(1), 2*x(2)])
%--------------------------------------------------------------------------
% AUTHOR: Sam Davenport
if nargin < 3
    fprime = NaN;
end
if nargin < 4
    tol = 10^(-8);
end
if nargin < 5
    max_iters = 100;
end

if size(initial_guess, 1) == 1
    initial_guess = initial_guess';
end
% warning('Using a different tolerance to the R Newton Raphson, may want to change this!!')

if ~isa(f, 'function_handle')
    error('f must be function handles. Eg @functionname')
end
fprimenan = 0;
if ~isa(fprime, 'function_handle')
    fprimenan = 1;
    if ~isnan(fprime)
        error('fprime must be NaN or a function handle')
    end
end

D = length(initial_guess);
try
    f(initial_guess);
    if ~fprimenan
        fprime(initial_guess);
    end
catch
    error('The input dimensions of the functions and the initial guess do not match')
end
if ~fprimenan && length(fprime(initial_guess)) ~= D %Tests that fprime(initial_guess) can be run while also testing the below
    error('The output number of dimensions of fprime should be the same as the dimension of the initial guess')
end

if isnan(sum(f(initial_guess)))
    error('f is not well defined')
end
if isnan(sum(fprime(initial_guess)))
    error('fprime is not well defined')
end

x_estimate = initial_guess;

iters = 0;
difference = sum(abs(f(x_estimate))) + 1; %plus 1 to ensure that it does actually go about the loop
% maxsizef = max(size(f));

if fprimenan
    while difference > tol
        gradmate = deriv_est(x_estimate, f);
        if sum(abs(gradmate(:))) < 100*eps
            warning('The derivative has reached a zero value')
            x_estimate = NaN;
            return
        end
        h = -pinv(gradmate)*f(x_estimate); %could change to: h_d = -f(x_estimate)/deriv_at_xest(d);
        x_estimate = x_estimate + h';
        
        difference = sum(abs(f(x_estimate))); %abs(f(new_x)) + abs(new_x - x_estimate)
        
        iters = iters + 1;
        if iters > max_iters
            error(strcat('The algorithm doesn''t convergence within ', num2str(max_iters), ' iterations'));
        end
    end
else
    while difference > tol
        gradmate = fprime(x_estimate);
        if sum(abs(gradmate(:))) < 100*eps
            warning('The derivative has reached a zero value')
            x_estimate = NaN;
            return
        end
        h = -pinv(gradmate)*f(x_estimate); %could change to: h_d = -f(x_estimate)/deriv_at_xest(d);
        x_estimate = x_estimate + h;
        
        difference = sum(abs(f(x_estimate))); %abs(f(new_x)) + abs(new_x - x_estimate)
        
        iters = iters + 1;
        if iters > max_iters
            error(strcat('The algorithm doesn''t convergence within ', num2str(max_iters), ' iterations'));
        end
    end
end

end