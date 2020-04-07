function [fderivplus, fderiv2] = derivest( fval, h, fplush, fminush )
% DERIVEST( fval, h, fplush, fminush )
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
%--------------------------------------------------------------------------
D = length(fplush);
fderivplus = zeros(1,D);
for I = 1:D
   fderivplus(I) = (fplush(I) - fval)/h;
end

if nargin > 3
    fderivminus = zeros(1,D);
    for I = 1:D
        fderivminus(I) = (fminush(I) - fval)/h;
    end
    fderiv = (1/2)*(fderivplus + fderivminus);
    fderiv2 = (fderivplus - fderivminus)/h;
else
    fderiv = fderivplus;
end


end

