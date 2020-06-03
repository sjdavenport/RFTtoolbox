%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%    This script tests the convfield_struct.m function
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% prepare workspace
clear all
close all

addpath( genpath( "/home/drtea/matlabToolboxes/RFTtoolbox/" ) )

% useful function to check equality of two arrays
sameArray = @(x,y) max(abs( x(:) - y(:) ));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Test section D=1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%% show that the Gaussian kernel does the same as using the naive method
%%%% without asuming seperable kernels
% Field parameters
FWHM = 10;
D = 1;
Nsubj = 120;
T = 100;
Dim = ones( [ 1 D ] ) * T;
siz =  ceil( 4* FWHM2sigma( FWHM ) );
mask = ones( [ Dim 1 ] );
sM = size(mask);

% resolution parameters
resAdd = 1;
dx = 1 / ( resAdd + 1 );
Dimhr = ( sM - 1 ) * resAdd + sM;

% generate data
Y = randn( [ Dim(1) Nsubj ]);
% increase the resolution of the raw data by introducing zeros
Y2 = zeros( [ Dimhr(1), Nsubj ] );
Y2( 1:( resAdd + 1 ):end, : ) = Y;

%%%% naive generation of convolution field
% grid for convolution kernel
x   = -siz:dx:siz;
% convolution kernel and derivatives to be used with convn
[ h, dxh ] = Gker( x, FWHM, 1 );     

% get the convolutional field (old)
smY_old  = convn( Y2, h', 'same' );
% get the derivative of the convolutional field
smYx_old = convn( Y2, dxh', 'same' );

%%%% convolution fields function
[ smY, xvals_vecs ] = convfield_struct( Y, FWHM, resAdd, D, 0, 1 );
DsmY = convfield_struct( Y, FWHM, resAdd, D, 1 );

%%%% awsome that still works!
[ sameArray( smY_old, smY ),...
sameArray( smYx_old, DsmY(:,:,:,1) ) ]

%%%% Plot slices of the image
close all
figure(1), clf,
subplot(1,2,1);
imagesc(smY_old), colorbar
title('naive')
subplot(1,2,2);
imagesc(smY), colorbar
title('convfield_struct')

% plot the different versions of computation of the fields. They seem to
% agree except for my naive use of fconv, so we can use it for the LKCestim
figure(2), clf,
subplot(1,2,1);
imagesc(smYx_old), colorbar
title('naive')
subplot(1,2,2);
imagesc(DsmY), colorbar
title('convfield_struct')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Test section D=2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%% show that the Gaussian kernel does the same as using the naive method
%%%% without asuming seperable kernels
% Field parameters
FWHM = 10;
D = 2;
Nsubj = 120;
T = 100;
Dim = ones( [ 1 D ] ) * T;
siz =  ceil( 4* FWHM2sigma( FWHM ) );
mask = ones( Dim );
sM = size(mask);

% resolution parameters
resAdd = 1;
dx = 1 / ( resAdd + 1 );
Dimhr = ( sM - 1 ) * resAdd + sM;

% generate data
Y = randn( [ Dim Nsubj ]);
% increase the resolution of the raw data by introducing zeros
Y2 = zeros( [ Dimhr, Nsubj ] );
Y2( 1:( resAdd + 1 ):end, 1:( resAdd + 1 ):end, : ) = Y;

%%%% naive generation of convolution field
% grid for convolution kernel
[x,y] = meshgrid( -siz:dx:siz, -siz:dx:siz );
xvals = [x(:), y(:)]';

% convolution kernels to be used ith convn
h   = reshape( GkerMV( xvals, FWHM ), size(x) );
dh  = GkerMVderiv( xvals, FWHM );
dxh = reshape( dh(2,:), size(x) );
dyh = reshape( dh(1,:), size(x) );

% get the convolutional field (old)
smY_old  = convn( Y2, h, 'same' );
% get the derivative of the convolutional field
smYx_old = convn( Y2, dxh, 'same' );
smYy_old = convn( Y2, dyh, 'same' );

%%%% convolution fields function
[ smY, xvals_vecs ] = convfield_struct( Y, FWHM, resAdd, D, 0, 1 );
DsmY = convfield_struct( Y, FWHM, resAdd, D, 1 );

%%%% awsome that still works!
[ sameArray( smY_old, smY ),...
sameArray( smYx_old, DsmY(:,:,:,1) ),...
sameArray( smYy_old, DsmY(:,:,:,2) ) ]

%%%% Plot slices of the image
close all
figure(1), clf,
subplot(1,2,1);
imagesc(smY_old(:,:,1)), colorbar
title('naive')
subplot(1,2,2);
imagesc(smY(:,:,1)), colorbar
title('convfield_struct')

% plot the different versions of computation of the fields. They seem to
% agree except for my naive use of fconv, so we can use it for the LKCestim
figure(2), clf,
subplot(1,2,1);
imagesc(smYx_old(:,:,1)), colorbar
title('naive')
subplot(1,2,2);
imagesc(DsmY(:,:,1,1)), colorbar
title('convfield_struct')

figure(3), clf,
subplot(1,2,1);
imagesc(smYy_old(:,:,1)), colorbar
title('naive')
subplot(1,2,2);
imagesc(DsmY(:,:,1,2)), colorbar
title('convfield_struct')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%    Test section D=2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%% show that the Gaussian kernel does the same as using the naive method
%%%% without asuming seperable kernels
% Field parameters
FWHM = 3;
D = 3;
Nsubj = 120;
T = 40;
Dim = ones( [ 1 D ] ) * T;
siz =  ceil( 4* FWHM2sigma( FWHM ) );
mask = ones( Dim );
sM = size(mask);

% resolution parameters
resAdd = 1;
dx = 1 / ( resAdd + 1 );
Dimhr = ( sM - 1 ) * resAdd + sM;

% generate data
Y = randn( [ Dim Nsubj ]);
% increase the resolution of the raw data by introducing zeros
Y2 = zeros( [ Dimhr, Nsubj ] );
Y2( 1:( resAdd + 1 ):end, 1:( resAdd + 1 ):end, 1:( resAdd + 1 ):end, : ) = Y;

%%%% naive generation of convolution field
% grid for convolution kernel
[x, y, z] = meshgrid( -siz:dx:siz, -siz:dx:siz, -siz:dx:siz );
xvals = [x(:), y(:), z(:)]';

% convolution kernels to be used ith convn
h   = reshape( GkerMV( xvals, FWHM ), size(x) );
dh  = GkerMVderiv( xvals, FWHM );
dxh = reshape( dh(2,:), size(x) );
dyh = reshape( dh(1,:), size(x) );
dzh = reshape( dh(3,:), size(x) );

% get the convolutional field (old)
smY_old  = convn( Y2, h, 'same' );
% get the derivative of the convolutional field
smYx_old = convn( Y2, dxh, 'same' );
smYy_old = convn( Y2, dyh, 'same' );
smYz_old = convn( Y2, dzh, 'same' );

%%%% convolution fields function
[ smY, xvals_vecs ] = convfield_struct( Y, FWHM, resAdd, D );
DsmY = convfield_struct( Y, FWHM, resAdd, D, 1 );

%%%% not as perfect as before, but should be good enough. probably errors
%%%% are adding up somehow???
[ sameArray( smY_old, smY ),...
  sameArray( smYx_old, DsmY(:,:,:,:,1) ),...
  sameArray( smYy_old, DsmY(:,:,:,:,2) ),...
  sameArray( smYz_old, DsmY(:,:,:,:,3) ) ]

%%%% Plot slices of the image
close all
figure(1), clf,
subplot(1,2,1);
imagesc(smY_old(:,:,20,1)), colorbar
title('naive')
subplot(1,2,2);
imagesc(smY(:,:,20,1)), colorbar
title('convfield_struct')

% plot the different versions of computation of the fields. They seem to
% agree except for my naive use of fconv, so we can use it for the LKCestim
figure(2), clf,
subplot(1,2,1);
imagesc(smYx_old(:,:,20,1)), colorbar
title('naive')
subplot(1,2,2);
imagesc(DsmY(:,:,20,1,1)), colorbar
title('convfield_struct')

figure(3), clf,
subplot(1,2,1);
imagesc(smYy_old(:,:,20,1)), colorbar
title('naive')
subplot(1,2,2);
imagesc(DsmY(:,:,20,1,2)), colorbar
title('convfield_struct')

figure(4), clf,
subplot(1,2,1);
imagesc(smYz_old(:,:,20,1)), colorbar
title('naive')
subplot(1,2,2);
imagesc(DsmY(:,:,20,1,3)), colorbar
title('convfield_struct')