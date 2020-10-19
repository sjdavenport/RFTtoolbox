function [ Lambda, FWHM_est ] = stat_Lambda( field, dfield )
% stat_Lambda( cfield, dcfield ) estimates the Reimmeinian metric under the 
% assumption of stationarity
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  field   an object of class Field representing observations of a field,
%          fiberD = 1 and fibersize > 1.
%  dfield  an object of class Field representing observations of the 
%          derivatives of a field.      
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow, Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
D = field.D;

%% Main function
%--------------------------------------------------------------------------
% Construct VoxManifold object by providing Riemannian metric
g = Riemmetric_est( field, dfield );
g = Mask( g );
G = reshape( g.field, [ prod( g.masksize ), D, D ] );

% Calculate Lambda by averaging over the mask
Lambda = squeeze( mean( G( g.mask(:), :, : ) ) );

% Apply the scaling factor
Lambda = Lambda*(field.fibersize-3)/(field.fibersize-2);
FWHM_est = sqrt(4*log(2)./diag(Lambda));

return
