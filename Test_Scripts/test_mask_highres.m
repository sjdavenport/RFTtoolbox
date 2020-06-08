%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%    This script tests the mask_highres.m function
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all

addpath( genpath( "/home/drtea/matlabToolboxes/RFTtoolbox/" ) )

% always show educational plots of the mask
plots = 1;


%% %-----------------------------------------------------------------------
%    Test section D = 1
%--------------------------------------------------------------------------
% resolution added
resAdd  = 1;
% enlargement of the original mask in high res voxels
enlarge = 0;

% create a mask and show it
Sig = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
mask = logical( Sig > 0.02 );
mask = mask(50,:)';
plot( mask ),
title( 'mask on orignal resolution' )
clear Sig

% create the high res mask and show the educational figures
mask_hr = mask_highres( mask, resAdd, enlarge, plots );


%% %-----------------------------------------------------------------------
%    Test section D = 2
%--------------------------------------------------------------------------
close all
% resolution added
resAdd  = 1;
% enlargement of the original mask in high res voxels
enlarge = ceil(resAdd/2);

% create a mask and show it
Sig = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
mask = logical( Sig > 0.02 & Sig < 1.1 );
figure, clf,
imagesc( mask ), colorbar,
title( 'mask on orignal resolution' )
clear Sig

% non enlarged domain
[ mask_hr, weights ] = mask_highres( mask, resAdd, 0, plots );
% enlarged domain
mask_hr = mask_highres( mask, resAdd, enlarge, plots );


%% % simpler example
close all
% resolution added
resAdd  = 3;
% enlargement of the original mask in high res voxels
enlarge = ceil(resAdd/2);

% create a mask and show it
mask = logical( ones(50,50) );
mask(:,50) = 0;
figure, clf,
imagesc( mask ), colorbar,
title( 'mask on orignal resolution' )

% non enlarged domain
mask_hr = mask_highres( mask, resAdd, 0, plots );
% enlarged domain
mask_hr = mask_highres( mask, resAdd, enlarge, plots );

%% %-----------------------------------------------------------------------
%    Test section D = 3
%--------------------------------------------------------------------------
close all
% resolution added
resAdd  = 3;
% enlargement of the original mask in high res voxels
enlarge = 1;

% create a mask and show it
siz = 13;
dx = 0.25;
[x,y,z] = meshgrid( -siz:dx:siz, -siz:dx:siz, -siz:dx:siz );
xvals = [x(:), y(:), z(:)]';
h   = reshape( GkerMV( xvals, 5 ), size(x) );
mask    = logical( h > 0.002 );
imagesc( mask(:,:,53) ),
colorbar,
title( 'mask on orignal resolution' )
clear h

% enlarged domain
mask_hr = mask_highres( mask, resAdd, enlarge, 0 );
% non enlarged domain
mask_hr = mask_highres( mask, resAdd, 0, 0 );