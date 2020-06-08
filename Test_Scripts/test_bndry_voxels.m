%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the bndry_voxels function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Test section D = 1
%--------------------------------------------------------------------------
close all
% create a mask and show it
Sig = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
mask = logical( Sig > 0.02 );
mask = mask(50,:)';
plot( mask ),
title( 'mask' )
clear Sig

% note in 1D there is only the 'full' option
bdry = bdry_voxels( mask, "full" );
figure(1), clf,
plot( mask + bdry )
title( 'mask + mask of boundary'  )

%% Test section D = 2
%--------------------------------------------------------------------------
%%%%%% Example demonstrating behaviour for all ones mask
% 2D examples - full
mask = ones(10,10);
bdry_voxels( logical(mask), 'full' )

% 2D examples - sides
bdry_voxels( logical(mask), 'x' )
bdry_voxels( logical(mask), 'y' )

%%%%%% simple box example
close all
% create a simple mask and show it
mask = zeros( [ 15 15 ] );
mask( 4:12, 5:10 ) = 1;
mask = logical( mask );
figure(1), clf,
imagesc( mask(:,:) ),
colorbar
title( 'mask' )

%%% plot the different options of boundary
% fixed y option
bdry = bdry_voxels( mask, "x" );
figure(2), clf,
imagesc( mask + bdry ), colorbar
title("boundary for fixed y directions")

% fixed x option
bdry = bdry_voxels( mask, "y" );
figure(3), clf,
imagesc( mask + bdry ), colorbar
title("boundary for fixed x directions")

% fixed "full" option
bdry = bdry_voxels( mask, "full" );
figure(4), clf,
imagesc( mask + bdry ), colorbar
title("all boundary points")

%% %%%% complicated mask example
close all
% generate mask
Sig = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
mask = logical( Sig > 0.02 & Sig < 1.1 );
figure(1), clf,
imagesc( mask ), colorbar,
title("mask")
clear Sig

%%% plot the different options of boundary
% fixed y option
bdry = bdry_voxels( mask, "x" );
figure(2), clf,
imagesc( mask + bdry ), colorbar
title("boundary for fixed y directions")

% fixed x option
bdry = bdry_voxels( mask, "y" );
figure(3), clf,
imagesc( mask + bdry ), colorbar
title("boundary for fixed x directions")

% fixed "full" option
bdry = bdry_voxels( mask, "full" );
figure(4), clf,
imagesc( mask + bdry ), colorbar
title("all boundary points")

%% Test section D = 3
%--------------------------------------------------------------------------
close all
% create a mask and show it
mask = zeros( [ 5 5 5 ] );
mask( 2:4, 2:4, 2:4 ) = 1;
figure(1), clf,
imagesc( mask(:,:,2) ), colorbar
title("slice of mask")

% get boundary voxels lying in z value planes or better having a not to the
% mask connected face pointing into the z direction
bdry = bdry_voxels( logical( mask ), "xy" )

% same as before for y
bdry = bdry_voxels( logical( mask ), "xz" )

% same as before for x
bdry = bdry_voxels( logical( mask ), "yz" )
