function FWHM = sigma2FWHM( sigma )
% FWHM = sigma2FWHM( sigma ) kernel bandwidths into FWHM.
%--------------------------------------------------------------------------
% ARGUMENTS
%   sigma: numeric vector/array of bandwidths for a smoother
%--------------------------------------------------------------------------
% OUTPUT
%   FWHM: the full width half maximum of the istropic kernel with bandwidth
%         sigma
%--------------------------------------------------------------------------
% EXAMPLES
% sigma = sigma2FWHM( 3 )
% FWHM = FWHM2sigma(sigma)
%--------------------------------------------------------------------------
% AUTHORS: Fabian Telschow and Samuel Davenport
%--------------------------------------------------------------------------
FWHM = sigma * sqrt( 8*log(2) );
end