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
FWHM   = 20;
pad    = ceil( 4*FWHM2sigma( FWHM ) );
% Method = "numerical";
method = "analytical";

%% Example with recangular mask
% Generate mask
mask = pad_vals( true( [ T, 1 ] ), pad, 0 );

% LKC from continuous theory
theoryL = LKC_isogauss_theory( FWHM, T  );

% Generate test data
lat_data = randn( [ T+2*pad nsubj ] );

% Closest approximation of the continuous field uses thresholding
% mask_lat = 0, since otherwise there are boundary effects
mask_lat = 0;
LKC1 = LKC_conv_est( lat_data, mask, FWHM, 1, mask_lat );
LKC3 = LKC_conv_est( lat_data, mask, FWHM, 3, mask_lat );
LKC5 = LKC_conv_est( lat_data, mask, FWHM, 5, mask_lat );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
[ theoryL; LKC1.hatL; LKC3.hatL; LKC5.hatL ]'

% Masking the data shows a boundary effect, validating departure
% from the theoretical value
LKC_masked = LKC_conv_est( lat_data, mask, FWHM, 1, 1 );
[ theoryL; LKC_masked.hatL ]'

%% %% D = 2 
%% % Parameters for the field
T      = 49;
nsubj  = 100;
FWHM   = sigma2FWHM(5);
pad = ceil( 4*FWHM2sigma( FWHM ) );

%% Example with recangular mask
% Get mask
mask = pad_vals( true( [ T T ] ), pad );

% Get LKC for the theoretical field
theoryL = LKC_isogauss_theory( FWHM, [ T T ] );

% Generate test data
lat_data = randn( [ T+2*pad T+2*pad nsubj ] );

% Closest approximation of the continuous field uses thresholding
% mask_lat = 0, since otherwise there are boundary effects
mask_lat = 0;
LKC1 = LKC_conv_est( lat_data, mask, FWHM, 1, mask_lat );
LKC3 = LKC_conv_est( lat_data, mask, FWHM, 3, mask_lat );
LKC5 = LKC_conv_est( lat_data, mask, FWHM, 5, mask_lat );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
[ theoryL; LKC1.hatL; LKC3.hatL; LKC5.hatL ]'

%% Example with complicated mask
Sig = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
mask = logical( Sig > 0.02 & Sig < 1.1 );
figure(1), clf,
imagesc( mask ), colorbar,
title("mask")
clear Sig

% Get LKC for the theoretical field
theoryL = LKC_wncfield_theory( mask, FWHM, 3, mask_lat );

% Generate test data
lat_data = randn( [ size(mask) nsubj ] );

% Closest approximation of the continuous field uses thresholding
% mask_lat = 0, since otherwise there are boundary effects
mask_lat = 0;
LKC1 = LKC_conv_est( lat_data, mask, FWHM, 1, mask_lat );
LKC3 = LKC_conv_est( lat_data, mask, FWHM, 3, mask_lat );
LKC5 = LKC_conv_est( lat_data, mask, FWHM, 5, mask_lat );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
[ theoryL; LKC1.hatL; LKC3.hatL; LKC5.hatL ]'

% Masking the data gives different LKCs
theoryL_masked = LKC_wncfield_theory( mask, FWHM, 3, 1 );
LKC_masked = LKC_conv_est( lat_data, mask, FWHM, 1, 1 );
[ theoryL; theoryL_masked; LKC_masked.hatL ]'

%% %% D = 3 
% Parameters for the field
T      = 20;
nsubj  = 10;
FWHM   = sigma2FWHM(1.5);
pad    = ceil( 4*FWHM2sigma( FWHM ) );

%% Rectangular domain example
% Generate rectangular mask with a padded zero collar 
mask = pad_vals( ones( [ T T T] ), pad );

% Get theoretical LKC
theoryL = LKC_wncfield_theory( mask, FWHM, 3, 0 );

% Generate test data
lat_data = randn( [ T+2*pad T+2*pad T+2*pad nsubj ] );

% Closest approximation of the continuous field uses thresholding
% mask_lat = 0, since otherwise there are boundary effects
mask_lat = 0;
LKC1 = LKC_conv_est( lat_data, mask, FWHM, 1, mask_lat );
LKC3 = LKC_conv_est( lat_data, mask, FWHM, 3, mask_lat );
LKC5 = LKC_conv_est( lat_data, mask, FWHM, 5, mask_lat );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
[ theoryL; LKC1.hatL; LKC3.hatL; LKC5.hatL ]'