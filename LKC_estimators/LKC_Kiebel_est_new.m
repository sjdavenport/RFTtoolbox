function [ L, Lambda, FWHM ] = LKC_Kiebel_est_new( field, k, num_deriv )
% LKC_kiebel_est( field, k, num_deriv ) estimates the LKCs assuming isotropy
% according to the article [1].
%
% [1] Kiebel, Stefan J., et al. "Robust smoothness estimation in
% statistical parametric maps using standardized residuals from the
% general linear model." Neuroimage 10.6 (1999): 756-766.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  field   an object of class Field representing observations of a field,
%          fiberD = 1 and fibersize > 1.
% Optional
%  k      integer such that df = field.fieldsize - k. It is the total
%         degree of freedoms of the input Field object, which in general
%         should be residuals, e.g., from a linear model, compare Kiebel
%         1999 eq. (10).
%  num_deriv  a string. Choices are "symmetric" or "onesided" for symmetric
%             or onesided numeric derivatives.
%--------------------------------------------------------------------------
% OUTPUT
%  L       an 1 x field.D vector containing the LKCs L_1,...,L_field.D
%  Lambda  a D x D matrix of the estimated Riemannian metric. Note this is
%          assumed constant across the domain, since Kiebel assumes
%          stationarity.
%  FWHM    a 1 x D vector of estimates of FWHM based on the Lambda matrix 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow, Samuel Davenport
%--------------------------------------------------------------------------

%% Check mandatory input
%--------------------------------------------------------------------------

%% Main function
%--------------------------------------------------------------------------
% Estimate the Lambda matrix from the grid
Lambda =  Lambda_lat_stat_est( field, num_deriv, k );

% Get the FHWM estimate
FWHM =  sqrt( 4*log(2) ./ diag(Lambda) );
voxmfd  = VoxManifold( constfield( Lambda, field.mask, field.xvals ) );

% Obtain the LKC estimate
L = LKC_est( voxmfd );
end
