function [ FWHM, sigma ] = Gker2FWHM( kernel )
%  Gker2FWHM( kernel ) obtains the FWHM given an isotropic Gaussian kernel
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  kernel       an object of class SepKernel
%--------------------------------------------------------------------------
% OUTPUT
% FWHM      the FWHM used to smooth in each direction
% sigma     the parameter of the kernel used to smooth in each direction
%--------------------------------------------------------------------------
% EXAMPLES
% D = 2; FWHM = 3;
% kernel = SepKernel( D, FWHM );
% Gker2FWHM(kernel)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Main Function Loop
%--------------------------------------------------------------------------
Ker0 = kernel.kernel{1}(0);
sigma = 1/(Ker0*sqrt(2*pi));
FWHM = sigma2FWHM(sigma);

end

