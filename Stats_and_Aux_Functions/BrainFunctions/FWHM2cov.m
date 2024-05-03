function Cx = FWHM2cov( x, FWHM )
% FWHM2COV( FWHM, D ) generates the covariance function
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  FWHM
%  D
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% Cx = FWHM2cov( 0, 3 )
% Cx = FWHM2cov( 1, 3 )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
D = length(x);

%%  Main Function Loop
%--------------------------------------------------------------------------
Lambda = FWHM2Lambda( FWHM, D );
Cx = exp((-1/4)*x'*(2*Lambda)*x);

end

