%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%    This script tests the bdry_voxels.m function
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% prepare workspace
clear all
close all

addpath( genpath( "/home/drtea/matlabToolboxes/RFTtoolbox/" ) )

%% %%%% 1D
% create a mask and show it
Sig = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
mask = logical( Sig > 0.02 );
mask = mask(50,:)';
figure(1), clf,
plot( mask )
clear Sig

bdry = bdry_voxels( mask, "full" );
figure(1), clf,
plot( mask + bdry )

%% %%%% 2D
close all
% create a simple mask and show it
mask = zeros( [ 5 5 ] );
mask( 2:4, 2:4 ) = 1;
mask = logical( mask );
figure(1), clf,
imagesc( mask(:,:) ), colorbar

bdry = bdry_voxels( mask, "x" );
figure(2), clf,
imagesc( mask + bdry ), colorbar
title("boundary for fixed y directions")

bdry = bdry_voxels( mask, "y" );
figure(3), clf,
imagesc( mask + bdry ), colorbar
title("boundary for fixed x directions")

bdry = bdry_voxels( mask, "full" );
figure(4), clf,
imagesc( mask + bdry ), colorbar
title("all boundary points")

% create a complicated mask and show it
Sig = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
mask = logical( Sig > 0.02 & Sig < 1.1 );
figure(5), clf,
imagesc( mask ), colorbar
clear Sig

bdry = bdry_voxels( mask, "x" );
figure(6), clf,
imagesc( mask + bdry ), colorbar
title("boundary for fixed y directions")

bdry = bdry_voxels( mask, "y" );
figure(7), clf,
imagesc( mask + bdry ), colorbar
title("boundary for fixed x directions")

bdry = bdry_voxels( mask, "full" );
figure(8), clf,
imagesc( mask + bdry ), colorbar
title("all boundary points")

%% %%%% 3D
% create a mask and show it
mask = zeros( [ 5 5 5 ] );
mask( 2:4, 2:4, 2:4 ) = 1;
figure(1), clf,
imagesc( mask(:,:,2) ), colorbar


bdry = bdry_voxels( logical( mask ), "xy" );

figure(2), clf,
imagesc( mask + bdry ), colorbar

bdry = bdry_voxels( mask, "y" );
figure(3), clf,
imagesc( mask + bdry ), colorbar

bdry = bdry_voxels( mask, "full" );
figure(4), clf,
imagesc( mask + bdry ), colorbar
