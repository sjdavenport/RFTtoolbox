function Lambda = FWHM2Lambda( FWHM, D )
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
% Lambda = FWHM2Lambda( 3, 2 )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

Lambda = 4*log(2)*eye(D)/FWHM^2;

end

