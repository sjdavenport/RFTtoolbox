function x_estimate = findpeak( initial_guess, f, fprime, fprime2, max_or_min, tol, max_iters)
% FINDPEAK( initial_guess, f, fprime, fprime2, max_or_min, tol, max_iters)
% finds the peak (either a maximum or minimum) of a random field by
% searching from an initial estimate.
%--------------------------------------------------------------------------
% ARGUMENTS
% initial_guess     an initial guess of the peak location
% f                 a function handle representing the function
% fprime            a function handle representing the function's derivative
%                   note that if fprime is NaN it is estimated using f.
% fprime2           a function handle representing the function's 2nd derivative
%                   note that if fprime2 is NaN it is estimated using fprime
% max_or_min        0/1 specifies whether to search for a local maximum or
%                   minimum. Default is 1: local maximum.
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
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
if nargin < 4
    fprime2 = NaN;
end
if nargin < 5
    max_or_min = 1;
end
if nargin < 6
    tol = 10^(-8);
end
if nargin < 7
    max_iters = 100;
end
if max_or_min == 0
    max_or_min = -1;
end
htol = 10^(-3);

if size(initial_guess, 1) == 1
    initial_guess = initial_guess';
end
% warning('Using a different tolerance to the R Newton Raphson, may want to change this!!')

if ~isa(f, 'function_handle')
    error('fprime must be function handles. Eg @functionname')
end
if ~isa(fprime, 'function_handle')
    error('fprime must be function handles. Eg @functionname')
end
fprime2nan = 0;
if ~isa(fprime2, 'function_handle')
    fprime2nan = 1;
    if ~isnan(fprime2)
        error('fprime2 must be NaN or a function handle')
    end
end

D = length(initial_guess);
try
    f(initial_guess);
    fprime(initial_guess);
    if ~fprime2nan
        fprime2(initial_guess);
    end
catch
    error('The input dimensions of the functions and the initial guess do not match')
end
if ~fprime2nan && length(fprime2(initial_guess)) ~= D %Tests that fprime(initial_guess) can be run while also testing the below
    error('The output number of dimensions of fprime should be the same as the dimension of the initial guess')
end

if isnan(sum(f(initial_guess)))
    error('f is not well defined')
end
if isnan(sum(fprime(initial_guess)))
    error('fprime is not well defined')
end
if isnan(sum(fprime2(initial_guess)))
    error('fprime2 is not well defined')
end

x_estimate = initial_guess;

iters = 0;
difference = sum(abs(fprime(x_estimate))) + 1; %plus 1 to ensure that it does actually go about the loop
% maxsizef = max(size(f));

if fprime2nan
    while difference > tol
        gradmate = deriv_est(x_estimate, fprime);
        if sum(abs(gradmate(:))) < 100*eps
            warning('The derivative has reached a zero value')
            x_estimate = NaN;
            return
        end
        h = -pinv(gradmate)*fprime(x_estimate); %could change to: h_d = -f(x_estimate)/deriv_at_xest(d);
        x_estimate = x_estimate + h';
        
        difference = sum(abs(fprime(x_estimate))); %abs(f(new_x)) + abs(new_x - x_estimate)
        
        iters = iters + 1;
        if iters > max_iters
            x_estimate = NaN*ones(1, D);
        end
    end
else
    while difference > tol
        gradmate = fprime2(x_estimate);
        if sum(abs(gradmate(:))) < 100*eps
            warning('The derivative has reached a zero value')
            x_estimate = NaN*ones(1, D);
            return
        end
        
        fprintf('f\n')
        f(x_estimate)
        fprintf('fprime\n')
        fprime(x_estimate)
        fprintf('eigenvalues of fprime2\n')
        eig(gradmate)
        fprintf('h\n')
        h = -pinv(gradmate)*fprime(x_estimate);
       
        if f(x_estimate + h) > f(x_estimate - h)
            x_estimate = x_estimate + max_or_min*h;
        else
            x_estimate = x_estimate - max_or_min*h;
        end
        f(x_estimate)
        
        [V,D] = eig(gradmate);
        edirections = find(max_or_min*diag(D) > 0);
        if all(h < htol) && ~isempty(edirections)
            x_estimate = x_estimate + 0.01*V(:,edirections(1));
            f(x_estimate)
        end
        
        difference = sum(abs(fprime(x_estimate))); %abs(f(new_x)) + abs(new_x - x_estimate)
        
        iters = iters + 1;
        if iters > max_iters
            x_estimate = NaN*ones(1, D);
        end
    end
end

end