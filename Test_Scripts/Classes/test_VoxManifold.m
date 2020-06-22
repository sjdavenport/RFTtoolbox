%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the VoxManifold class object and its functions
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%% Basic constructor
%% %% nargin = 1
%% % Input can be a complete fields object containing the Riemannian metric
Mfd = VoxManifold( wnfield( [ 3 4 5 ], [3, 3] ) )
%% % Input can be a vector defining the size of the mask.
% mask field is defaulted to be true( varargin{1} ) and xvals is
% defaulted by a sequence. Hence only a field need to be added.
Mfd = VoxManifold( [ 4 5 ] )
Mfd.g = wnfield( [ 4 5 ], [2 2] );
% mask or xvals can be changed, if needed. Note that masked property will
% change to 0, since the field is nonzero/NaN/-Inf outside the mask.
mask = true( [ 2, 3 ] )
mask = logical( pad_vals( mask ) )
Mfd.mask = mask

%% % Input can be a 1xD cell array containing vectors.
% mask field is defaulted to be of compatible dimension with the xvals.
% Hence only a field need to be added.
Mfd = VoxManifold({ 1:3, [ 0.5, 3, 6, 6.5 ], 15:19 })

%% % Input can be a logical mask
% xvals field is defaulted to be of compatible dimension with the xvals.
% Hence only a field need to be added.
mask = true( [ 2, 3 ] )
mask = logical( pad_vals( mask ) )
Mfd = VoxManifold( mask )

%% %% nargin = 2
%% % Input can be a mask and a xvals.
% xvals is defaulted, yet can be changed
mask = true( [ 4, 12 ] )
mask = logical( pad_vals( mask ) )

Mfd = VoxManifold( mask, { 1:6, (1:14)-5 } )

%% % Input can be a Riemannian metric g and resadd.
% xvals is defaulted, yet can be changed
mask = true( [ 4, 12 ] );
mask = logical( pad_vals( mask ) );
g    = EuclideanMetric( Field( mask ) );

Mfd = VoxManifold( g, 3 )

%% %% nargin = 3
%% % Input can specify a mask, xvals and resadd.
mask = true( [ 4, 12 ] )
mask = logical( pad_vals( mask ) )

Mfd = VoxManifold( mask, { 1:6, (1:14)-5 }, 3 )
%% % Input can be a Riemannian metric g, resadd and enlarge. 
% xvals is defaulted, yet can be changed
mask = true( [ 4, 12 ] );
mask = logical( pad_vals( mask ) );
g    = EuclideanMetric( Field( mask ) );

Mfd = VoxManifold( g, 3 )