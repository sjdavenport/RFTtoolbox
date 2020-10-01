function [ FWHM, Lambda ] = cov2FWHM( rho, D )
% cov2FWHM( rho ) calculates the FWHM from a value of the correlation
% between two adjacent voxels in a given dimension.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  rho      the correlation between two adjacent voxels
%  D        the number of dimensions
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% FWHM = 3;
% rho = FWHM2cov( 1, FWHM )
% cov2FWHM( rho, 1 )
% 
% cov2FWHM( rho^(1/100), 1)
% FWHM2cov( 1, 30 )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Main Function Loop
%--------------------------------------------------------------------------
for d = 1:D
    Lambda(d,d) = -2*log(rho);
end

FWHM = Lambda2FWHM(Lambda);

end

