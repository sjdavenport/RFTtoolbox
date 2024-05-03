function [Lambda, detLambda] = Lambda_theory( FWHM, D )
% LAMBDA_THEORY( FWHM, D )
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if length(FWHM) == 1
    FWHM = repmat(FWHM,1,D);
elseif length(FWHM) ~= D
    error('FWHM must have length 1 or D')
end

%%  Main Function Loop
%--------------------------------------------------------------------------
Gker_params = FWHM2sigma(FWHM);
Lambda      = diag(Gker_params.^(-2) / 2);
detLambda   = det(Lambda);
end

