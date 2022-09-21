function [ L, L0, nonstatInt ] = LKC_voxmfd_est( field, dfield, d2field, version, scale )
% LKC_voxmfd_est( field, dfield ) estimates the LKCs of the voxel
% manifold defined by the domain (aka its mask) of a convolution field and
% its induced Riemannian metric
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
%          - D = 3, logical of length 3. version(1) indicates whether L2
%                   should be estimated, version(2) whether the first 
%                   integral is used in L1 and version(3) whether also
%                   the second integral is used in L1. Default: [1 1 0];
%                   i.e., the stationary approximation of L1
%  scale    0/1 determines whether to scale by the correction factor or not
%--------------------------------------------------------------------------
% OUTPUT
% L   an 1 x field.D vector containing the LKCs L_1,...,L_field.D
% L0  a numeric containing the zeroth LKC, i.e. the Euler characteristic of
%     the domain.
%--------------------------------------------------------------------------
% EXAMPLES
% See test_LKC_voxmfd_est.m
%--------------------------------------------------------------------------
% Author: Fabian Telschow and Samuel Davenport
%--------------------------------------------------------------------------

%% Check and add optional input
%--------------------------------------------------------------------------

if ~exist( 'version', 'var' )
    if field.D < 3
        version = true;
    else
        version = logical( [ 1 1 0 ] );
    end
end

if ~exist('scale','var')
    scale = 0;
end

%% Main function
%--------------------------------------------------------------------------

% Construct VoxManifold object by providing Riemannian metric
Lambda = Riemmetric_est( field, dfield );

% Adjust by the scaling factor if specified
if scale == 1
    Lambda = Lambda*((field.fibersize-3)/(field.fibersize-2));
end
voxmfd = VoxManifold( Lambda );

if exist( 'd2field', 'var' ) && ~isempty(d2field)
    if ~isempty( d2field.field )
        voxmfd.Gamma = Christoffel_est( field, dfield, d2field );
    end
end

% Obtain the LKCs
[ L, L0, nonstatInt ] = LKC_est( voxmfd, version );

return