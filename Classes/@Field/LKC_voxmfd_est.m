function [ L, L0 ] = LKC_voxmfd_est( field, dfield, d2field, version )
% LKC_voxmfd_est( field, dfield ) estimates the LKCs of the voxel
% manifold defined by the domain (aka its mask) of a convolution field and
% its induced Riemannian metric.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  field   an object of class Field representing observations of a field,
%          fiberD = 1 and fibersize > 1.
%  dfield  an object of class Field representing observations of the 
%          derivatives of a field.      
%
% Optional
%  d2field  an object of class Field representing observations of the 
%           second derivatives of a field.
%  version a logical/ logical vector. Length depends on voxmfd.D
%          - D = 1, always true.
%          - D = 2, either 1 or 0. 1 if L1 should be estimated, 0 else.
%          - D = 3, logical of length 3. version(1), indicates whether L2
%          should be estimated, version(2) whether the first integral is
%          used in L1 and version(3) whether the second integral is used.
%
%--------------------------------------------------------------------------
% OUTPUT
% L   an 1 x field.D vector containing the LKCs L_1,...,L_field.D
% L0  a numeric containing the zeroth LKC, i.e. the Euler characteristic of
%     the domain.
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% Author: Fabian Telschow
%--------------------------------------------------------------------------

%% Check and add optional input
%--------------------------------------------------------------------------

if ~exist( 'version', 'var' )
    if field.D < 3
        version = true;
    else
        version = true( [ 1 3 ] );
    end
end

%% Main function
%--------------------------------------------------------------------------

% Construct VoxManifold object by providing Riemannian metric
voxmfd = VoxManifold( Riemmetric_est( field, dfield ) );

if exist( 'd2field', 'var' )
    if ~isempty( d2field.field )
        voxmfd.Gamma = Christoffel_est( field, dfield, d2field );
    end
end

% Obtain the LKCs
[ L, L0 ] = LKC_est( voxmfd, version );

return