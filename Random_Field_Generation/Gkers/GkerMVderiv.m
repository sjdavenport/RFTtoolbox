function deriv = GkerMVderiv( x, sigma2_or_FWHM, use_fwhm )
% GkerMVderiv( x, sigma2_or_FWHM, use_fwhm ) calculates the Gaussian Kernel given
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
% EXAMPLE
% val = GkerMV([1,2]', 3)
% 
% deriv = GkerMVderiv([1,2]', 3)
% h = 0.00001;
% valplushx = GkerMV([1+h,2]', 3);
% valplushx
% valplushy = GkerMV([1,2+h]', 3);
% (valplushx - val)/h
% (valplushy - val)/h
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
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
val = exp(-sum(x.^2,1)/(2*sigma2))/(sqrt(2*pi*sigma2)^D); 
deriv = -(Sigmainv*x).*repmat(val,D,1);
% deriv = -(Sigmainv*x).*val; Not used for 2016 compatibility

end

