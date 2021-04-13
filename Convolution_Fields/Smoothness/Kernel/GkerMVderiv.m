function deriv = GkerMVderiv( x, sigma2_or_FWHM, use_fwhm )
% GkerMVderiv( x, sigma2_or_FWHM, use_fwhm ) calculates the partial 
% derivatives of theGaussian Kernel given multivariate data and the FWHM.
%--------------------------------------------------------------------------
% ARGUMENTS
% x                 a D by nevals matrix where each column is a
%                   D-dimensional vector at which to evaluate the kernel.
% sigma2_or_FWHM    If FWHM, it is the FWHM in voxels.
% use_fwhm          0/1 specifies whether to use the FWHM or sigma2. Default
%                   is 1 i.e. to use the FWHM.
%--------------------------------------------------------------------------
% OUTPUT
% deriv         a D by nevals matrix whose jth column is the vector of
%               partial derivatives of the Gaussian Kernel with given FWHM 
%               evaluated at the jth column of x.
%--------------------------------------------------------------------------
% EXAMPLE
% % 2D
% val = GkerMV([1,2]', 3)
% deriv = GkerMVderiv([1,2]', 3)
% h = 0.00001;
% valplushx = GkerMV([1+h,2]', 3);
% valplushy = GkerMV([1,2+h]', 3);
% (valplushx - val)/h
% (valplushy - val)/h
%
% % 3D
% val = GkerMV([1,2,0.5]', 3)
% deriv = GkerMVderiv([1,2,0.5]', 3)
% h = 0.00001;
% valplushx = GkerMV([1+h,2,0.5]', 3);
% valplushy = GkerMV([1,2+h,0.5]', 3);
% valplushz = GkerMV([1,2,0.5+h]', 3);
% (valplushx - val)/h
% (valplushy - val)/h
% (valplushz - val)/h
%
% % 3D multiple points:
% GkerMVderiv([1,2,0.5;1,0.5,2]', 3)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
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

