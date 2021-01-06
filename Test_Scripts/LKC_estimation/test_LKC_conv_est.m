%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the LKC_conv_est.m function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% D = 1 
%% % Field parameters
T      = 100;
nsubj  = 120;
FWHM   = 2;
pad    = ceil( 4*FWHM2sigma( FWHM ) );

%% Example with recangular mask
% Generate mask
mask = pad_vals( true( [ T, 1 ] ), pad, false );

% Generate params object for convfields
resThy = 301;
params = ConvFieldParams( FWHM, resThy, ceil(resThy/2), false );

% LKC from continuous theory
theoryLt = LKC_isogauss_theory( FWHM, T  ); % This is for continuous fields
theoryL  = LKC_wncfield_theory( mask, params );

% Generate test data
lat_data = wnfield( mask, nsubj );

% Estimate across different resadd
params = ConvFieldParams( FWHM, 1, ceil(1/2), false );
LKC1   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( FWHM, 3, ceil(3/2), false );
LKC3   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( FWHM, 5, ceil(5/2), false );
LKC5   = LKC_latconv_est( lat_data, params );

% Values are stable accross different resadd increases.
% Note that resadd should be odd.Otherwise the underlying
% manifold changes
[ theoryL; LKC1; LKC3; LKC5 ]'

%% Example with recangular mask without enlargement
% Generate mask
mask = pad_vals( true( [ T, 1 ] ), pad, false );
FWHM = 3;

% Generate params object for convfields
resThy = 101
params = ConvFieldParams( FWHM, resThy, 0, false );

% LKC from continuous theory
theoryLt = LKC_isogauss_theory( FWHM, T  ); % This is for continuous fields
theoryL  = LKC_wncfield_theory( mask, params );

% Generate test data
lat_data = wnfield( mask, nsubj );

% Estimate across different resadd
params = ConvFieldParams( FWHM, 0, 0, false );
LKC0   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( FWHM, 1, 0, false );
LKC1   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( FWHM, 2, 0, false );
LKC2   = LKC_latconv_est( lat_data, params );

% Values are stable accross different resadd increases.
% Note that resadd should be odd.Otherwise the underlying
% manifold changes
[ theoryL; LKC0; LKC1; LKC2 ]'

%% D = 1 no pad vals example
% Generate mask
mask = true( [ T, 1 ] );

% Generate params object for convfields
resThy = 301
params = ConvFieldParams( FWHM, resThy, ceil(resThy/2), false );

% LKC from continuous theory
theoryLt = LKC_isogauss_theory( FWHM, T  );
theoryL  = LKC_wncfield_theory( mask, params );

% Generate test data
lat_data = wnfield( mask, nsubj );

% Estimate across different resadd
params = ConvFieldParams( FWHM, 1, ceil(1/2), false );
LKC1   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( FWHM, 3, ceil(3/2), false );
LKC3   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( FWHM, 5, ceil(5/2), false );
LKC5   = LKC_latconv_est( lat_data, params );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
[ theoryL; LKC1; LKC3; LKC5 ]'

%% D = 1 non-stationary sphere example
FWHM = 3;
% Dimension of domain
dim   = [ 100 1 ];

% Generate mask
pad         = ceil( 4 * FWHM2sigma(FWHM) );
mask        = true( dim );
mask(5:95) = false;
mask        = logical( pad_vals( mask, pad) );

% Mask the lattice data
lat_mask = false % switch to true for non-stationarity

% Generate params object for convfields
resThy = 101
params = ConvFieldParams( FWHM, resThy, 0, lat_mask );

% LKC from continuous theory
theoryLt = LKC_isogauss_theory( FWHM, T  );
theoryL  = LKC_wncfield_theory( mask, params );

theoryL

% Generate test data
lat_data = wnfield( mask, nsubj );

% Estimate across different resadd
params = ConvFieldParams( FWHM, 1, 0, lat_mask );
LKC1   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( FWHM, 3, 0, lat_mask );
LKC3   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( FWHM, 5, 0, lat_mask );
LKC5   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( FWHM, 15, 0, lat_mask );
LKC15   = LKC_latconv_est( lat_data, params );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
[ theoryL; LKC1; LKC3; LKC5; LKC15 ]'

%% %% D = 2 
%% % Parameters for the field
T      = 49;
nsubj  = 100;
FWHM   = sigma2FWHM(1.2);
pad = ceil( 4*FWHM2sigma( FWHM ) );

%% Example with recangular mask
% Get mask
mask = pad_vals( true( [ T T ] ), pad, false );

% Generate params object for convfields
params = ConvFieldParams( [ FWHM, FWHM ], 21, ceil(21/2), false );

% LKC from continuous theory
theoryLt = LKC_isogauss_theory( FWHM, [ T T ]  );
theoryL  = LKC_wncfield_theory( mask, params );

% Generate test data
lat_data = wnfield( mask, nsubj );

% Estimate across different resadd
params = ConvFieldParams( [ FWHM, FWHM ], 1, ceil(1/2), false );
LKC1   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( [ FWHM, FWHM ], 3, ceil(3/2), false );
LKC3   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( [ FWHM, FWHM ], 5, ceil(5/2), false );
LKC5   = LKC_latconv_est( lat_data, params );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
[ theoryL; LKC1; LKC3; LKC5 ]'

%% Example with recangular mask no enlargement
% Get mask
mask = pad_vals( true( [ T T ] ), pad, false );

% Generate params object for convfields
params = ConvFieldParams( [ FWHM, FWHM ], 21, 0, false );

% LKC from continuous theory
theoryLt = LKC_isogauss_theory( FWHM, [ T T ]  );
theoryL  = LKC_wncfield_theory( mask, params );

% Generate test data
lat_data = wnfield( mask, nsubj );

% Estimate across different resadd
params = ConvFieldParams( [ FWHM, FWHM ], 0, 0, false );
LKC0   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( [ FWHM, FWHM ], 1, 0, false );
LKC1   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( [ FWHM, FWHM ], 2, 0, false );
LKC2   = LKC_latconv_est( lat_data, params );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
[ theoryL; LKC0; LKC1; LKC2 ]'

%% Example with complicated mask
Sig = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
mask = logical( Sig > 0.02 & Sig < 1.1 );
figure(1), clf,
imagesc( mask ), colorbar,
title("mask")
clear Sig

% Generate params object for convfields
params = ConvFieldParams( [ FWHM, FWHM ], 21, ceil(21/2), false );

% LKC from continuous theory
theoryL  = LKC_wncfield_theory( mask, params );

% Generate test data
lat_data = wnfield( mask, nsubj );

% Estimate across different resadd
params = ConvFieldParams( [ FWHM, FWHM ], 1, ceil(1/2), false );
LKC1   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( [ FWHM, FWHM ], 3, ceil(3/2), false );
LKC3   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( [ FWHM, FWHM ], 5, ceil(5/2), false );
LKC5   = LKC_latconv_est( lat_data, params );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
[ theoryL; LKC1; LKC3; LKC5 ]'

%% %% D = 3 
% Parameters for the field
T      = 5;
nsubj  = 10;
FWHM   = sigma2FWHM(1.5);
pad    = ceil( 4*FWHM2sigma( FWHM ) );

%% Rectangular domain example
% Generate rectangular mask with a padded zero collar 
mask = pad_vals( ones( [ T T T ] ), [2 1 3], false );

% Get theoretical LKC
params = ConvFieldParams( [ FWHM, FWHM, FWHM ], 1, ceil(1/2), false );
theoryL  = LKC_wncfield_theory( mask, params );
% LKC from continuous theory
theoryLt = LKC_isogauss_theory( FWHM, [ T T T ]  );

% Generate test data
lat_data = wnfield( mask, nsubj );

% Estimate across different resadd
params = ConvFieldParams( [ FWHM, FWHM, FWHM ], 1, ceil(1/2), false );
LKC1   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( [ FWHM, FWHM, FWHM ], 3, ceil(3/2), false );
LKC3   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( [ FWHM, FWHM, FWHM ], 5, ceil(5/2), false );
LKC5   = LKC_latconv_est( lat_data, params );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
[ theoryL; LKC1; LKC3; LKC5 ]'

%% Testing 3D (there's a bug)
Dim = 2*[16,16,16]; mask = true(Dim); resadd = 1;
spfn = @(nsubj) normrnd(0,1,[Dim, nsubj]);
lat_data = spfn(10);
LKC_conv_est( lat_data, mask, FWHM, 1)

%% Pad vals vs no pad vals (comment/uncomment lines 148,149)
T      = 5;
nsubj  = 10;
FWHM   = sigma2FWHM(1.5);
pad    = ceil( 4*FWHM2sigma( FWHM ) );

% Generate rectangular mask with a padded zero collar 
% mask = pad_vals( ones( [ T T T] ), pad );
mask = true([T+2*pad T+2*pad T+2*pad]);

% Get theoretical LKC
theoryL = LKC_wncfield_theory( mask, FWHM, 3, 0 );
contL   = LKC_isogauss_theory( FWHM, [ T T T] );

% Generate test data
lat_data = randn( [ T+2*pad T+2*pad T+2*pad nsubj ] );

% Closest approximation of the continuous field uses thresholding
% mask_lat = 0, since otherwise there are boundary effects
mask_lat = 0;
LKC1 = LKC_conv_est( lat_data, mask, FWHM, 1 );