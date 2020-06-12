%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the bndry_voxels function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% D = 1 
%% Create a mask and show it
Sig = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
mask = logical( Sig > 0.02 );
mask = mask(50,:)';
plot( mask ),
title( 'mask' )
clear Sig

% Note that in 1D there is only the "full" option
bdry = bndry_voxels( mask', "full" );
figure(1), clf,
plot( mask + bdry )
title( 'mask + mask of boundary'  )

%% %% Test section D = 2
%% 2D examples - full - example demonstrating behaviour for all ones mask
mask = ones(10,10);
[bndry, weights] = bndry_voxels( logical(mask), "full" )

%% 2D examples - sides
[bndry, weights] = bndry_voxels( logical(mask), "x" )
[bndry, weights] = bndry_voxels( logical(mask), "y" )

%% Simple box example
% Create a simple mask and show it
mask = zeros( [ 15 15 ] );
mask( 4:12, 5:10 ) = 1;
mask = logical( mask );
figure(1), clf,
imagesc( mask(:,:) ),
colorbar
title( 'mask' )

%% % plot the different options of boundary
%% Fixed x option
bdry = bndry_voxels( mask, "x" );
figure(2), clf,
imagesc( mask + bdry ), colorbar
title("boundary for fixed x directions")

%% Fixed y option
bdry = bndry_voxels( mask, "y" );
figure(3), clf,
imagesc( mask + bdry ), colorbar
title("boundary for fixed y directions")

%% Fixed "full" option
bdry = bndry_voxels( mask, "full" );
figure(4), clf,
imagesc( mask + bdry ), colorbar
title("all boundary points")

%% Example demonstrating the bdry fix
% Create a simple mask and show it
mask = true( [ 10 10 ] );
mask( 1:3, 1:4 ) = 0;
mask = logical( mask );
figure(1), clf,
imagesc( mask(:,:) ),
colorbar
title( 'mask' )

%% % plot the different options of boundary
%% Fixed y option
bdry = bndry_voxels( mask, "y" );
figure(2), clf,
imagesc( mask + bdry ), colorbar
title("boundary for fixed y directions")

%% Fixed x option
bdry = bndry_voxels( mask, "x" );
figure(3), clf,
imagesc( mask + bdry ), colorbar
title("boundary for fixed x directions")

%% Fixed "full" option
bdry = bndry_voxels( mask, "full" );
figure(4), clf,
imagesc( mask + bdry ), colorbar
title("all boundary points")

%% Complicated mask example
% Generate mask
Sig = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
mask = logical( Sig > 0.02 & Sig < 1.1 );
figure(1), clf,
imagesc( mask ), colorbar,
title("mask")
clear Sig

%% % plot the different options of boundary
%% Fixed y option
bdry = bndry_voxels( mask, "y" );
figure(2), clf,
imagesc( mask + bdry ), colorbar
title("boundary for fixed y directions")

%% Fixed x option
bdry = bndry_voxels( mask, "x" );
figure(3), clf,
imagesc( mask + bdry ), colorbar
title("boundary for fixed x directions")

%% Fixed "full" option
bdry = bndry_voxels( mask, "full" );
figure(4), clf,
imagesc( mask + bdry ), colorbar
title("all boundary points")

%% Fixed x+y value
xbdry = bndry_voxels( mask, "x" );
ybdry = bndry_voxels( mask, "y" );
figure(3), clf,
imagesc( (ybdry & xbdry) + xbdry ), colorbar
title("boundary for fixed x directions")

%% %% Test section D = 3
%% % Simple usage example on a box
% Create a mask and show it
mask = zeros( [ 5 5 5 ] );
mask( 2:4, 2:4, 2:4 ) = 1;
figure(1), clf,
imagesc( mask(:,:,2) ), colorbar
title("slice of mask")

% Get all boundary voxels seperate by type
[ bdry, weights ] = bndry_voxels( logical( mask ) )

%% % Example demonstrating its use for high resolution masks
% Create a mask and show it
mask = true( [ 3 3 3 ] );
mask( 1, 1, 1 ) = 0;
visualize_bndry3D( mask, 1, ["x"], 40, [ -102 16 ] )
title("Resolution increased mask + x-edges")

resadd = 1;
mask_hr = mask_highres( mask, resadd, ceil(resadd/2) );

% Get all boundary voxels seperate by type
[ bdry_xy, weights_xy ] = bndry_voxels( logical( mask_hr ) )

% Visualize all edges
visualize_bndry3D( mask, 1, [ "x", "y", "z" ], 40, [ -102 16 ] )
title("All edge points from bndry_voxels")

% Visualize all edges
visualize_bndry3D( mask, 1, [ "xy", "yz", "xz" ], 40, [ -102 16 ] )
title("All faces points from bndry_voxels")
