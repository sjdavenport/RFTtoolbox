function FWHM = Lambda2FWHM( Lambda )
% FWHM2Lambda( FWHM, D ) computes the continuous theoretical Lambda matrix
% that results from smoothing white noisewith an isotropic Gaussian kernel 
% in D dimensions with the specified FWHM.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  FWHM     the FWHM of the isotropic Gaussian smoothing kernel
%  D        the number of dimensions
%--------------------------------------------------------------------------
% OUTPUT
% Lambda    the covariance matrix of the partial derivatives
%--------------------------------------------------------------------------
% EXAMPLES
% D = 2; FWHM = 3;
% Lambda = FWHM2Lambda( FWHM, D )
% FWHM = Lambda2FWHM(Lambda)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

D = size(Lambda,1);
FWHM = 1/sqrt((trace(Lambda)/D)/(4*log(2)));

end

