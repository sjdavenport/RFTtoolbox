function FWHM = sigma2FWHM( sigma )
% FWHM = sigma2FWHM( sigma )
% This function converts bandwidths into FWHM.
%--------------------------------------------------------------------------
% ARGUMENTS
%   sigma numeric vector/array of bandwidths for a smoother
%--------------------------------------------------------------------------
% OUTPUT
%   FWHM bandwidth converted to FWHM
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow.

FWHM = sigma * sqrt( 8*log(2) );
end