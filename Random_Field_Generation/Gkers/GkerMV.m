function [val, deriv, deriv2] = GkerMV( x, sigma2_or_FWHM, use_fwhm )
% GKERMV( x, sigma2_or_FWHM, use_fwhm ) calculates the Gaussian Kernel given
% multivariate data and the FWHM of the kernel.
%--------------------------------------------------------------------------
% ARGUMENTS
% x                 a D by nevals matrix where each column is a
%                   D-dimensional vector at which to evaluate the kernel.
% sigma2_or_FWHM    If FWHM, it is the FWHM in voxels.
% use_fwhm
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% [val, deriv, deriv2] = GkerMV([1.5,2], 3)
% 
% [val, deriv] = GkerMV([1,2]', 3)
% h = 0.00001;
% valplushx = GkerMV([1+h,2]', 3); h = 0.00001;
% valplushx
% valplushy = GkerMV([1,2+h]', 3); h = 0.00001;
% (valplushx - val)/h
% (valplushy - val)/h
% [~, derivplushx, deriv2] = GkerMV([1+h,2]', 3);
% [~, derivplushy, deriv2] = GkerMV([1,2+h]', 3);
% deriv2 
% (derivplushx - deriv)/h
% (derivplushy - deriv)/h
%
% [val, deriv, deriv2] = GkerMV([1,2]', 3)
% [val, deriv, deriv2] = GkerMV([1,3;2,4], 3)
%--------------------------------------------------------------------------
% AUTHOR: Samuel J. Davenport
% Don't change this! Too much depends on it! LOL
if nargin < 3
    use_fwhm = 1;
end

if use_fwhm
    sigma2 = FWHM2sigma(sigma2_or_FWHM)^2;
else
    sigma2 = sigma2_or_FWHM;
end

D = size(x, 1);

Sigmainv = (1/sigma2)*eye(D);
val = exp(-sum(x.^2,1)/(2*sigma2))/(sqrt(2*pi*sigma2)^D); %Only this it is Multivariate!
temp_deriv = -Sigmainv*x; %Without the kernel constant.

deriv2 = zeros(D^2, size(x,2));
for I = 1:D
    for J = 1:I
        deriv2( D*(I-1) + J, : ) = (temp_deriv(I,:).*temp_deriv(J,:) - Sigmainv(I,J)).*val;
        deriv2( D*(J-1) + I, : ) = deriv2( D*(I-1) + J, : );
    end
end 
deriv = temp_deriv.*val;

%For the moment just have an isotropic kernel coded. Need to generalize!

end

