function [ tcf_deriv, tcf_deriv2 ] = MVderiv( field, tval )
% MVderiv(field, h) calculates the derivative and second derivative of the
% field at given tvals.
%--------------------------------------------------------------------------
% ARGUMENTS
% tval      the t values (an ndim=D by nvalues matrix) at which to evaluate 
%           the derivatives.
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
if nargin < 3
    h = 0.00001;
end
if h > 0.00001
    warning('Small h is bad due to numerical imprecision')
end
D = size(tval, 1);

if D == 1
    tcf_deriv = @(tval) (field(tval+h) - field(tval))/h;
    tcf_deriv2 = @(tval) (tcf_deriv(tval+h) - tcf_deriv(tval))/h;
elseif D == 2
    tcf_deriv = @(x) [(field(x+h*[1,0]') - field(x))/h, (field(x+h*[0,1]') - field(x))/h]';
    tcf_deriv2 = @(x) [(tcf_deriv(x+h*[1,0]') - tcf_deriv(x))/h, (tcf_deriv(x+h*[0,1]') - tcf_deriv(x))/h];
elseif D == 3
    tcf_deriv = @(x) [(field(x+h*[1,0,0]') - field(x))/h, (field(x+h*[0,1,0]') - field(x))/h, (field(x+h*[0,0,1]') - field(x))/h]';
    tcf_deriv2 = @(x) [(tcf_deriv(x+h*[1,0,0]') - tcf_deriv(x))/h, (tcf_deriv(x+h*[0,1,0]') - tcf_deriv(x))/h, (tcf_deriv(x+h*[0,0,1]') - tcf_deriv(x))/h];
end

end

