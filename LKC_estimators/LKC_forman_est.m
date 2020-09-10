function [ L_forman, L_kiebel, FWHM_f, FWHM_k ] = LKC_forman_est( field, df )
% LKC_forman_est( field ) estimates the LKCs assuming isotropy according
% to the article CITE.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  field   an object of class Field representing observations of a field,
%          fiberD = 1 and fibersize > 1.
% Optional
%  df      integer representing the degrees of freedom of the field. This
%          is needed, if the residuals are 
%--------------------------------------------------------------------------
% OUTPUT
%  L   an 1 x field.D vector containing the LKCs L_1,...,L_field.D
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR:  Samuel Davenport, Fabian Telschow
%--------------------------------------------------------------------------

%% Check mandatory input
%--------------------------------------------------------------------------

%% Main function
%--------------------------------------------------------------------------
% Riemannian metric estimated using Forman's formula
[ FWHM_f, FWHM_k, Lambda_k, ~] =  est_smooth_field( field, df );

Lambda_f = 4*log(2) * diag( 1 ./ FWHM_f.^2 );

% Obtain the LKCs from Forman FWHM
voxmfd  = VoxManifold( constfield( Lambda_f, field.mask, field.xvals ) );
L_forman = LKC_est( voxmfd );

% Obtain the LKCs from Forman FWHM
voxmfd  = VoxManifold( constfield( Lambda_k, field.mask, field.xvals ) );
L_kiebel = LKC_est( voxmfd );