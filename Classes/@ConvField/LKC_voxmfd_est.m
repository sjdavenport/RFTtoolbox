function [ L, L0 ] = LKC_voxmfd_est( cfield, dcfield )
% LKC_voxmfd_est( obj1, obj2, masked ) estimates the LKCs of the voxel
% manifold defined by a convolution field.
% The convolution fields should be generated using convfield.m.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  cfield   an object of class ConvField of derivtype 0 and fiberD = 1 and
%           fibersize > 1.
%  dcfield  an object of class ConvField of derivtype 1 and fiberD = 1 and
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

% Construct VoxManifold object by providing Riemannian metric and resadd
% and enlarge
voxmfd = VoxManifold( Lambda_est( cfield, dcfield ),...
                      cfield.resadd, cfield.enlarge );               

% Obtain the LKCs
[ L, L0 ] = LKC_est( voxmfd );

return