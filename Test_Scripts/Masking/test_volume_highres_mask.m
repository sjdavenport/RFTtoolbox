%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests volumes of parts of the boundary on different
%%%    resolutions
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function to compute the volume
get_volume = @( weights, resadd, D ) sum( weights(:) ) * 1/(resadd+1)^D;

%%
%%% 1D
D = 1
% create a mask and show it
Sig  = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
mask = logical( Sig > 0.02 );
mask = mask(50,:)';
clear Sig

% weights for resolution added
[~, weights1, old1] = mask_highres( mask, 1, ceil(1/2) );
[~, weights3, old3] = mask_highres( mask, 3, ceil(3/2) );
[~, weights5, old5] = mask_highres( mask, 5, ceil(5/2) );

[[ get_volume( mask, 0, D ), get_volume( weights1, 1, D ),...
   get_volume( weights3, 3, D ), get_volume( weights5, 5, D )];
 [ get_volume( mask, 0, D ), get_volume( old1, 1, D ),...
   get_volume( old3, 3, D ), get_volume( old5, 5, D ) ]]

%% 1D with boundary
mask = true(1,10);
[~, weights1, old1] = mask_highres( mask, 1, ceil(1/2) );
[~, weights3, old3] = mask_highres( mask, 3, ceil(3/2) );
[~, weights5, old5] = mask_highres( mask, 5, ceil(5/2) );
[[ get_volume( mask, 0, D ), get_volume( weights1, 1, D ),...
   get_volume( weights3, 3, D ), get_volume( weights5, 5, D )];
 [ get_volume( mask, 0, D ), get_volume( old1, 1, D ),...
   get_volume( old3, 3, D ), get_volume( old5, 5, D ) ]]

%%
%%% 2D
D = 2
% create a mask and show it
Sig  = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
mask = logical( Sig > 0.02 & Sig < 1.1 );
imagesc( mask )
clear Sig

% weights for resolution added
[mask1, weights1, old1] = mask_highres( mask, 1, ceil(1/2) );
[mask2, weights3, old3] = mask_highres( mask, 3, ceil(3/2) );
[mask3, weights5, old5] = mask_highres( mask, 5, ceil(5/2) );

[[ get_volume( mask, 0, D ), get_volume( weights1, 1, D ),...
   get_volume( weights3, 3, D ), get_volume( weights5, 5, D )];
 [ get_volume( mask, 0, D ), get_volume( old1, 1, D ),...
   get_volume( old3, 3, D ), get_volume( old5, 5, D ) ]]

% Check volume of boundary in diffent resolutions
dmask  = bndry_voxels( mask,  "x" ) +  bndry_voxels( mask,  "y" );
dmask1 = bndry_voxels( mask1, "x" ) +  bndry_voxels( mask1, "y" );
dmask3 = bndry_voxels( mask2, "x" ) +  bndry_voxels( mask2, "y" );
dmask5 = bndry_voxels( mask3, "x" ) +  bndry_voxels( mask3, "y" );

[ get_volume( dmask, 0, D-1 ), get_volume( dmask1, 1, D-1 ),...
  get_volume( dmask3, 3, D-1 ), get_volume( dmask5, 5, D-1 )]

%% Bug, note comments in getweights
D = 2
% Here is something wrong
mask = true(3,3);
mask(1,2) = 0;

% weights for resolution added
[mask1, weights1, old1] = mask_highres( mask, 1, ceil(1/2) );
[mask2, weights3, old3] = mask_highres( mask, 3, ceil(3/2) );
[mask3, weights5, old5] = mask_highres( mask, 5, ceil(5/2) );

[[ get_volume( mask, 0, D ), get_volume( weights1, 1, D ),...
   get_volume( weights3, 3, D ), get_volume( weights5, 5, D )];
 [ get_volume( mask, 0, D ), get_volume( old1, 1, D ),...
   get_volume( old3, 3, D ), get_volume( old5, 5, D ) ]]

% Your example (Here not, so the boundary of the image is treated wrongly)
mask = false(5,5);
mask(2:4,2:4) = 1;
mask(2,3) = 0;

% weights for resolution added
[mask1, weights1, old1] = mask_highres( mask, 1, ceil(1/2) );
[mask2, weights3, old3] = mask_highres( mask, 3, ceil(3/2) );
[mask3, weights5, old5] = mask_highres( mask, 5, ceil(5/2) );

[[ get_volume( mask, 0, D ), get_volume( weights1, 1, D ),...
   get_volume( weights3, 3, D ), get_volume( weights5, 5, D )];
 [ get_volume( mask, 0, D ), get_volume( old1, 1, D ),...
   get_volume( old3, 3, D ), get_volume( old5, 5, D ) ]]


%%
%%% 3D
D  = 3;

% create a mask and show it
siz = 8;
dx  = 0.5;
[x,y,z] = meshgrid( -siz:dx:siz, -siz:dx:siz, -siz:dx:siz );
xvals = [x(:), y(:), z(:)]';
h     = reshape( GkerMV( xvals, 5 ), size(x) );
mask  = logical( h > 0.002 );
clear h

% weights for resolution added
[~, weights1, old1] = mask_highres( mask, 1, ceil(1/2) );
[~, weights3, old3] = mask_highres( mask, 3, ceil(3/2) );
[~, weights5, old5] = mask_highres( mask, 5, ceil(5/2) );

[[ get_volume( mask, 0, D ), get_volume( weights1, 1, D ),...
   get_volume( weights3, 3, D ), get_volume( weights5, 5, D )];
 [ get_volume( mask, 0, D ), get_volume( old1, 1, D ),...
   get_volume( old3, 3, D ), get_volume( old5, 5, D ) ]]

% Check volume of boundary in diffent resolutions
dmask  = bndry_voxels( mask,  "x" ) +  bndry_voxels( mask,  "y" );
dmask1 = bndry_voxels( mask1, "x" ) +  bndry_voxels( mask1, "y" );
dmask3 = bndry_voxels( mask2, "x" ) +  bndry_voxels( mask2, "y" );
dmask5 = bndry_voxels( mask3, "x" ) +  bndry_voxels( mask3, "y" );

[ get_volume( dmask, 0, D-1 ), get_volume( dmask1, 1, D-1 ),...
  get_volume( dmask3, 3, D-1 ), get_volume( dmask5, 5, D-1 )]


%%
D = 3
mask = true( 3, 3, 3 );
mask(1,1,1) = 0;

% weights for resolution added
[~, weights1, old1] = mask_highres( mask, 1, ceil( 1/2 ) );
[~, weights3, old3] = mask_highres( mask, 3, ceil( 3/2 ) );
[~, weights5, old5] = mask_highres( mask, 5, ceil( 5/2 ) );

[[ get_volume( mask, 0, D ), get_volume( weights1, 1, D ),...
   get_volume( weights3, 3, D ), get_volume( weights5, 5, D )];
 [ get_volume( mask, 0, D ), get_volume( old1, 1, D ),...
   get_volume( old3, 3, D ), get_volume( old5, 5, D ) ]]

%% Bug, note comments in getweights
% Here is something wrong
mask = true(3,3,3);
mask(1,1,2) = 0;

% weights for resolution added
[~, weights1, old1] = mask_highres( mask, 1, ceil(1/2) );
[~, weights3, old3] = mask_highres( mask, 3, ceil(3/2) );
[~, weights5, old5] = mask_highres( mask, 5, ceil(5/2) );

[[ get_volume( mask, 0, D ), get_volume( weights1, 1, D ),...
   get_volume( weights3, 3, D ), get_volume( weights5, 5, D )];
 [ get_volume( mask, 0, D ), get_volume( old1, 1, D ),...
   get_volume( old3, 3, D ), get_volume( old5, 5, D ) ]]

% Your example (Here not, so the boundary of the image is treated wrongly)
mask = false(5,5,5);
mask(2:4,2:4,2:4) = 1;
mask(2,2,3) = 0;

% weights for resolution added
[~, weights1, old1] = mask_highres( mask, 1, ceil(1/2) );
[~, weights3, old3] = mask_highres( mask, 3, ceil(3/2) );
[~, weights5, old5] = mask_highres( mask, 5, ceil(5/2) );

[[ get_volume( mask, 0, D ), get_volume( weights1, 1, D ),...
   get_volume( weights3, 3, D ), get_volume( weights5, 5, D )];
 [ get_volume( mask, 0, D ), get_volume( old1, 1, D ),...
   get_volume( old3, 3, D ), get_volume( old5, 5, D ) ]]

%% %% Boundary volume
D  = 3;
%% % Spherical object
% Create a mask and show it
siz = 3;
dx  = 0.5;
[x,y,z] = meshgrid( -siz:dx:siz, -siz:dx:siz, -siz:dx:siz );
xvals = [x(:), y(:), z(:)]';
h     = reshape( GkerMV( xvals, 5 ), size(x) );
mask  = logical( h > 0.003 );
imagesc( mask(:,:,7) )
clear h

%% Get resolution increased masks
bdry_vols = zeros([ 1 3 ]);
k = 0;
for resadd = [ 1 3 5 ]
    k = k+1;
    dx = 1 / ( resadd + 1 );
    mask_hr = mask_highres( mask, resadd, ceil(resadd/2) );
    [ bdry_xy, weights_xy ] = bndry_voxels( mask_hr, "xy" );
    [ bdry_xz, weights_xz ] = bndry_voxels( mask_hr, "xz" );
    [ bdry_yz, weights_yz ] = bndry_voxels( mask_hr, "yz" );

    bdry_vols(k) = ( sum( weights_xy(:) ) + ...
                  sum( weights_xz(:) ) + ...
                  sum( weights_yz(:) ) ) * dx^2;
end
bdry_vols

%% % Box object
% Create a mask and show it
mask  = true([ 10 5 10 ]);
clear h

%% Get resolution increased masks
bdry_vols = zeros([ 1 3 ]);
k = 0;
for resadd = [ 1 3 5 ]
    k = k+1;
    dx = 1 / ( resadd + 1 );
    mask_hr = mask_highres( mask, resadd, ceil(resadd/2) );
    [ bdry_xy, weights_xy ] = bndry_voxels( mask_hr, "xy" );
    [ bdry_xz, weights_xz ] = bndry_voxels( mask_hr, "xz" );
    [ bdry_yz, weights_yz ] = bndry_voxels( mask_hr, "yz" );

    bdry_vols(k) = ( sum( weights_xy(:) ) + ...
                  sum( weights_xz(:) ) + ...
                  sum( weights_yz(:) ) ) * dx^2;
end
bdry_vols

%% % Box object with corner missing
% Create a mask and show it
mask  = true([ 10 5 10 ]);
mask(1:3, 1:2, 1:3) = 0;
clear h

%% Get resolution increased masks
bdry_vols = zeros([ 1 3 ]);
bdry_vols_xy = zeros([ 1 3 ]);
bdry_vols_yz = zeros([ 1 3 ]);
bdry_vols_xz = zeros([ 1 3 ]);
k = 0;
for resadd = [ 1 3 5 ]
    k = k+1;
    dx = 1 / ( resadd + 1 );
    mask_hr = mask_highres( mask, resadd, ceil(resadd/2) );
    [ bdry_xy, weights_xy ] = bndry_voxels( mask_hr, "xy" );
    [ bdry_xz, weights_xz ] = bndry_voxels( mask_hr, "xz" );
    [ bdry_yz, weights_yz ] = bndry_voxels( mask_hr, "yz" );

    bdry_vols_xy(k) = sum( weights_xy(:) ) * dx^2;
    bdry_vols_xz(k) = sum( weights_xz(:) ) * dx^2;
    bdry_vols_yz(k) = sum( weights_yz(:) ) * dx^2;
end
bdry_vols_xy
bdry_vols_xz
bdry_vols_yz
bdry_vols_xy(k) + bdry_vols_yz + bdry_vols_xz;
