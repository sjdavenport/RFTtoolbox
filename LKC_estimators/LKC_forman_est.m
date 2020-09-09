function [ L, Lambda, FWHM ] = LKC_forman_est( field, df )
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
% Get the important variables from the Field object
Y    = field.field;
mask = zero2nan( field.mask );

% check that method is implemented for dimension D
if( D > 4 )
    error('The method is currently only implemented for field domains of dimension D<5.')
end

%% Main function
%--------------------------------------------------------------------------
% Riemannian metric estimated using Forman's formula
[ FWHM, ~, Lambda, ~] = est_smooth( Y, mask );

% Generate the voxel manifold
voxmfd  = VoxManifold( constfield( Lambda, field.mask, field.xvals ) );

% Obtain the LKCs
L = LKC_est( voxmfd );