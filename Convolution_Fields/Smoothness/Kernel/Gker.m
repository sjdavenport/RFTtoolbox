function [val, deriv, deriv2] = Gker( x, sigma2_or_FWHM, use_fwhm )
% GKER( x, sigma2_or_FWHM, use_fwhm ) calculates the Gaussian Kernel given
% data and the variance: sigma2 or FWHM.
%--------------------------------------------------------------------------
% ARGUMENTS
% x
% sigma2_or_FWHM    If FWHM, it is the FWHM in voxels.
% use_fwhm          Using the FWHM is the default
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
% AUTHOR: Samuel Davenport
if nargin < 3
    use_fwhm = 1;
end

if use_fwhm
    sigma2 = FWHM2sigma(sigma2_or_FWHM)^2;
else
    sigma2 = sigma2_or_FWHM;
end

val = exp(-x.^2/(2*sigma2))/sqrt(2*pi*sigma2);
deriv = (-x/sigma2).*exp(-x.^2/(2*sigma2))/sqrt(2*pi*sigma2); 
% dxh      = -x .* h / nu( 1 )^2;
deriv2 = (-1/sigma2 + x.^2/sigma2^2).*exp(-x.^2/(2*sigma2))/sqrt(2*pi*sigma2);

end

