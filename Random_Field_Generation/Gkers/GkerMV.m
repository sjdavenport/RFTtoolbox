function [val, deriv, deriv2] = GkerMV( x, sigma2_or_FWHM, use_fwhm )
% GKER( x, sigma2_or_FWHM, use_fwhm ) calculates the Gaussian Kernel given
% multivariate data and the FWHM of the kernel.
%--------------------------------------------------------------------------
% ARGUMENTS
% x
% sigma2_or_FWHM    If FWHM, it is the FWHM in voxels.
% use_fwhm
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% Gker([1.5,2]', 3)
% 
% Gker([1.5,2], 3)
% Gker([1.5,2; 0, 1], 3)
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

deriv = (-x/sigma2).*exp(-x.^2/(2*sigma2))/sqrt(2*pi*sigma2);
deriv2 = (-1/sigma2 + x.^2/sigma2^2).*exp(-x.^2/(2*sigma2))/sqrt(2*pi*sigma2);

D = size(x, 2);
val = exp(-sum(x.^2,2)/(2*sigma2))/(sqrt(2*pi*sigma2)^D);
%For the moment just have an isotropic kernel coded. Need to generalize!

end

