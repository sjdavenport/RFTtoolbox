function sigma = FWHM2sigma( FWHM )
% FWHM2sigma( FWHM ) converts the FWHM to std deviation
%--------------------------------------------------------------------------
% ARGUMENTS
% FWHM:     the FWHM of a smoothing kernel
%--------------------------------------------------------------------------
% OUTPUT
% sigma:    the value in terms of the field
%--------------------------------------------------------------------------
% EXAMPLES
% FWHM = FWHM2sigma( 5 )
% sigma = sigma2FWHM( FWHM )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
sigma = FWHM/sqrt(8*log(2));
end

