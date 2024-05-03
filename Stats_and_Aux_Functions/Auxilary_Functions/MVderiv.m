function [ f_deriv, f_deriv2 ] = MVderiv( f, D, h )
% MVderiv(field, h) calculates the derivative and second derivative of the
% field at given tvals.
%--------------------------------------------------------------------------
% ARGUMENTS
% f      a function handle that takes a t-value in R^D as input and
%        returns a scalar value.
% D      an integer representing the dimension of the input space of f.
% h      default is h = 0.00001
%--------------------------------------------------------------------------
% OUTPUT
% f_deriv     a function that takes t in R^D to \nabla f(t)
% f_deriv2    a function that takes t in R^D to \nabla^2 f(t)
%--------------------------------------------------------------------------
% EXAMPLES
% % 1D
% f = @(x) sin(x);
% derivfn = MVderiv(f, 1)
% derivfn(0)
% 
% f = @(x) cos(x);
% derivfn = MVderiv(f, 1)
% derivfn(0)
% 
% % 2D
% f = @(x) (x(1)^2 + x(2)^2);
% [derivfn, deriv2fn] = MVderiv(f, 2);
% derivfn([1,1]')
% deriv2fn([1,1]') 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if ~exist('h', 'var')
    h = 0.00001;
end
if h < 0.00001
    warning('Small h is bad due to numerical imprecision')
end

if D == 1
    f_deriv = @(tval) (f(tval+h) - f(tval))/h;
    f_deriv2 = @(tval) (f_deriv(tval+h) - f_deriv(tval))/h;
elseif D == 2
    f_deriv = @(x) [(f(x+h*[1,0]') - f(x))/h, (f(x+h*[0,1]') - f(x))/h]';
    f_deriv2 = @(x) [(f_deriv(x+h*[1,0]') - f_deriv(x))/h, (f_deriv(x+h*[0,1]') - f_deriv(x))/h];
elseif D == 3
    f_deriv = @(x) [(f(x+h*[1,0,0]') - f(x))/h, (f(x+h*[0,1,0]') - f(x))/h, (f(x+h*[0,0,1]') - f(x))/h]';
    f_deriv2 = @(x) [(f_deriv(x+h*[1,0,0]') - f_deriv(x))/h, (f_deriv(x+h*[0,1,0]') - f_deriv(x))/h, (f_deriv(x+h*[0,0,1]') - f_deriv(x))/h];
end

end

