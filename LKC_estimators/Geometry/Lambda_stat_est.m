function [ Lambda, FWHM_est ] = Lambda_stat_est( field, dfield, scale )
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
% FWHM = 6; resadd = 1; nsubj = 100; nvox = 100;
% params = ConvFieldParams(FWHM, resadd);
% trunc = ceil(4*FWHM2sigma(FWHM));
% noise = wfield(nvox + 2*trunc, nsubj);
% smooth_f = convfield(noise, params);
% start = floor((trunc/(nvox+2*trunc))*smooth_f.masksize(1));
% stat_f = smooth_f(start:smooth_f.masksize(1)-start);
% smooth_f_deriv = convfield(noise, params, 1);
% stat_f_derivs = smooth_f_deriv(start:smooth_f.masksize(1) -start);
% stat_Lambda( stat_f, stat_f_derivs )
% FWHM2Lambda(FWHM, 1)
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

% Adjust by the scaling factor, where does it come from?
if scale == 1
    Lambda = Lambda*(g.fibersize-3)/(g.fibersize-2);
end
FWHM_est = sqrt(4*log(2)./diag(Lambda));

return
