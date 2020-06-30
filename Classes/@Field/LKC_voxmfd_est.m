function [ L, L0 ] = LKC_voxmfd_est( cfield, dcfield, d2cfield, version )
% LKC_voxmfd_est( cfield, dcfield ) estimates the LKCs of the voxel
% manifold defined by the domain (aka its mask) of a convolution field and
% its induced Riemannian metric.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  cfield   an object of class ConvField of derivtype 0, fiberD = 1 and
%           fibersize > 1.
%  dcfield  an object of class ConvField of derivtype 1, fiberD = 1 and
%           fibersize > 1.
%--------------------------------------------------------------------------
% OUTPUT
% L   an 1 x cfield.D vector containing the LKCs L_1,...,L_cfield.D
% L0  a numeric containing the zeroth LKC, i.e. the Euler characteristic of
%     the domain.
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% Author: Fabian Telschow
%--------------------------------------------------------------------------

%% Main function
%--------------------------------------------------------------------------

% Construct VoxManifold object by providing Riemannian metric, resadd
% and enlarge
voxmfd = VoxManifold( Riemmetric_est( cfield, dcfield ),...
                      cfield.resadd, cfield.enlarge );
voxmfd.Gamma = Christoffel_est( cfield, dcfield, d2cfield );

% Obtain the LKCs
[ L, L0 ] = LKC_est( voxmfd, version );

return