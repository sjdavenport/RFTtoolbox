function [ L, Lambda, FWHM ] = LKC_Forman_est_new( field, k, num_deriv, pool )
%  LKC_Forman_est( field, k, num_deriv ) estimates the LKCs assuming isotropy according
% to the article [1].
%
% [1] Forman, Steven D., et al. "Improved assessment of significant
% activation in functional magnetic resonance imaging (fMRI): Use of a
% cluster‐size threshold." Magnetic Resonance in medicine 33.5 (1995): 636-647.
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
if ~exist('pool', 'var')
   pool = false; 
end

%% Main function
%--------------------------------------------------------------------------
% Estimate the Lambda matrix from the grid
Lambda =  Lambda_lat_stat_est( field, num_deriv, k, pool );

% Forman estimators
sigma_est = sqrt( -1 ./ 4 ./ log( 1 - diag(Lambda) / 2 ) );
FWHM      = sigma2FWHM(sigma_est);

% Get the Riemannian metric estimator
Lambda = 4 * log(2) * diag( 1 ./ FWHM.^2 );

% Obtain the LKCs from Forman Riemannian metric
voxmfd  = VoxManifold( constfield( Lambda, field.mask, field.xvals ) );
L = LKC_est( voxmfd );
%--------------------------------------------------------------------------



