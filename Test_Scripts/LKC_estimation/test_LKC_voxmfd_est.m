%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the LKC_voxmfd_est.m function
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
FWHM   = 10;
pad    = ceil( 4*FWHM2sigma( FWHM ) );

Mboot = 1e3;

%% Example with recangular mask
% Generate mask
mask = logical( pad_vals( true( [ T, 1 ] ), pad, 0 ) );

% LKC from continuous theory
theoryL = LKC_isogauss_theory( FWHM, T  );

% Generate test data
lat_data = wnfield( mask, nsubj );

% approximate continuous field
params = ConvFieldParams( FWHM, 1 );
params.lat_masked = false;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC1 = LKC_voxmfd_est( cfield, dcfield );

params = ConvFieldParams( FWHM, 3 );
params.lat_masked = false;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC3 = LKC_voxmfd_est( cfield, dcfield );

params = ConvFieldParams( FWHM, 5 );
params.lat_masked = false;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC5 = LKC_voxmfd_est( cfield, dcfield );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
[ theoryL; LKC1; LKC3; LKC5 ]'

%% D = 1 no pad vals example
% Calculate mask
mask = true( [ T, 1 ] );

% LKC from continuous theory
theoryL = LKC_isogauss_theory( FWHM, T  );

% Generate test data
lat_data = wnfield( mask, nsubj );

% approximate continuous field
params = ConvFieldParams( FWHM, 1 );
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC1 = LKC_voxmfd_est( cfield, dcfield );

params = ConvFieldParams( FWHM, 3 );
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC3 = LKC_voxmfd_est( cfield, dcfield );

params = ConvFieldParams( FWHM, 5 );
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC5 = LKC_voxmfd_est( cfield, dcfield );

% Values differ slightly from theory because of boundary effect.
[ theoryL; LKC1; LKC3; LKC5 ]'

%% %% D = 2 
%% % Parameters for the field
D      = 2
T      = 49;
nsubj  = 100;
FWHM   = sigma2FWHM(5);
pad = ceil( 4*FWHM2sigma( FWHM ) );

Mboot = 1e3;

%% Example with recangular mask
% Get mask
mask = logical( pad_vals( true( [ T T ] ), pad ) );

% Get LKC for the theoretical field
theoryL = LKC_isogauss_theory( FWHM, [ T T ] );

% Generate test data
lat_data = wnfield( mask, nsubj );

% approximate continuous field
params = ConvFieldParams( FWHM * ones([1 D]), 1 );
params.lat_masked = false;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC1 = LKC_voxmfd_est( cfield, dcfield );

params = ConvFieldParams( FWHM * ones([1 D]), 3 );
params.lat_masked = false;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC3 = LKC_voxmfd_est( cfield, dcfield );

params = ConvFieldParams( FWHM * ones([1 D]), 5 );
params.lat_masked = false;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC5 = LKC_voxmfd_est( cfield, dcfield );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
[ theoryL; LKC1; LKC3; LKC5 ]'

%% Example with complicated mask
Sig = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
mask = logical( Sig > 0.02 & Sig < 1.1 );
figure(1), clf,
imagesc( mask ), colorbar,
title("mask")
clear Sig

lat_masked = false;

% Get LKC for the theoretical field
params_thy = ConvFieldParams( FWHM * ones([1 D]), 7 );
params_thy.lat_masked = lat_masked; 
theoryL = LKC_wncfield_theory( mask, params_thy );

% Generate test data
lat_data = wnfield( mask, nsubj );

% approximate continuous field
params = ConvFieldParams( FWHM * ones([1 D]), 1 );
params.lat_masked = lat_masked;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC1 = LKC_voxmfd_est( cfield, dcfield );

params = ConvFieldParams( FWHM * ones([1 D]), 3 );
params.lat_masked = lat_masked;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC3 = LKC_voxmfd_est( cfield, dcfield );

params = ConvFieldParams( FWHM * ones([1 D]), 5 );
params.lat_masked = lat_masked;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC5 = LKC_voxmfd_est( cfield, dcfield );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
[ theoryL; LKC1; LKC3; LKC5 ]'

%% %% D = 3 
% Parameters for the field
D      = 3;
T      = 5;
nsubj  = 30;
FWHM   = 3;
pad    = 0;%ceil( 4*FWHM2sigma( FWHM ) );
lat_masked = true;
dim  = [ T T T];
dimp = dim + 2 * pad;

Mboot = 3e2; % You can increase it, yet it will take longer

%% Rectangular domain example
% Generate rectangular mask with a padded zero collar 
mask = logical( pad_vals( ones( dim ), pad ) );

% Get theoretical LKC
params_thy = ConvFieldParams( FWHM * ones([1 D]), 7 );
params_thy.lat_masked = lat_masked;
theoryL = LKC_wncfield_theory( mask, params_thy );
contL = LKC_isogauss_theory( FWHM, dim );

% Generate test data
lat_data = wnfield( mask, nsubj );

% approximate continuous field
params = ConvFieldParams( FWHM * ones([1 D]), 1 );
params.lat_masked = lat_masked;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC1 = LKC_voxmfd_est( cfield, dcfield );
LKC1_HPE = LKC_HP_est( cfield, Mboot, 1 );

params = ConvFieldParams( FWHM * ones([1 D]), 3 );
params.lat_masked = lat_masked;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC3 = LKC_voxmfd_est( cfield, dcfield );
LKC3_HPE = LKC_HP_est( cfield, Mboot, 1 );

params = ConvFieldParams( FWHM * ones([1 D]), 5 );
params.lat_masked = lat_masked;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC5 = LKC_voxmfd_est( cfield, dcfield );
LKC5_HPE = LKC_HP_est( cfield, Mboot, 1 );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
[ contL; theoryL; LKC1; LKC3; LKC5 ]'
[ contL; theoryL; LKC1_HPE.hatL'; LKC3_HPE.hatL'; LKC5_HPE.hatL' ]'

%% Sphere domain example
% Generate sphere mask with a padded zero collar 
mask = logical( pad_vals( ones( dim ), pad ) );
mask = bndry_voxels( mask, 'full' );

% Get theoretical LKC
params_thy = ConvFieldParams( FWHM * ones([1 D]), 7 );
params_thy.lat_masked = lat_masked;
theoryL = LKC_wncfield_theory( mask, params_thy );

% Generate test data
lat_data = wnfield( mask, nsubj );

% approximate continuous field
params   = ConvFieldParams( FWHM * ones([1 D]), 1 );
params.lat_masked = lat_masked;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC1     = LKC_voxmfd_est( cfield, dcfield );
LKC1_HPE = LKC_HP_est( cfield, Mboot, 1 );

params = ConvFieldParams( FWHM * ones([1 D]), 3 );
params.lat_masked = lat_masked;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC3 = LKC_voxmfd_est( cfield, dcfield );
LKC3_HPE = LKC_HP_est( cfield, Mboot, 1 );

params = ConvFieldParams( FWHM * ones([1 D]), 5 );
params.lat_masked = lat_masked;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC5 = LKC_voxmfd_est( cfield, dcfield );
LKC5_HPE = LKC_HP_est( cfield, Mboot, 1 );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
ConvE = [ theoryL; LKC1; LKC3; LKC5 ]'
HPE = [ theoryL; LKC1_HPE.hatL'; LKC3_HPE.hatL'; LKC5_HPE.hatL' ]'

%% Sphere domain example (highly non stationary data, L1 estimated using
% only "locally stationary" integral)
% Generate sphere mask with a padded zero collar 
mask = logical( pad_vals( ones( dim ), pad ) );
mask = bndry_voxels( mask, 'full' );

% Get theoretical LKC
params_thy.lat_masked = lat_masked;
theoryL = LKC_wncfield_theory( mask, params_thy );

% Generate test data
cut   = 6;
shift = 4;
voxmap  = 1:prod( dimp );
tmp2    = mod( voxmap, cut);
voxmap2 = [ voxmap( tmp2 <= cut-shift ), voxmap( tmp2 > cut-shift ) ];

lat_data = cnfield( mask, 1, voxmap2, 2, nsubj );

% approximate continuous field
params   = ConvFieldParams( FWHM * ones([1 D]), 1 );
params.lat_masked = lat_masked;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC1     = LKC_voxmfd_est( cfield, dcfield );
LKC1_HPE = LKC_HP_est( cfield, Mboot, 1 );

params = ConvFieldParams( FWHM * ones([1 D]), 3 );
params.lat_masked = lat_masked;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC3 = LKC_voxmfd_est( cfield, dcfield );
LKC3_HPE = LKC_HP_est( cfield, Mboot, 1 );

params = ConvFieldParams( FWHM * ones([1 D]), 5 );
params.lat_masked = lat_masked;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
LKC5 = LKC_voxmfd_est( cfield, dcfield );
LKC5_HPE = LKC_HP_est( cfield, Mboot, 1 );

% Here we see that the second integral somehow must be wrong
ConvE = [ theoryL; LKC1; LKC3; LKC5 ]'
HPE = [ theoryL; LKC1_HPE.hatL'; LKC3_HPE.hatL'; LKC5_HPE.hatL' ]'

%% Sphere domain example (highly non stationary data, L1 estimated using
% both integrals)
% Generate sphere mask with a padded zero collar 
mask = logical( pad_vals( ones( dim ), pad ) );
mask = bndry_voxels( mask, 'full' );

% Get theoretical LKC
params_thy = ConvFieldParams( FWHM * ones([1 D]), 7 );
theoryL = LKC_wncfield_theory( mask, params_thy );

% Generate test data
cut   = 6;
shift = 4;
voxmap  = 1:prod( dimp );
tmp2    = mod( voxmap, cut);
voxmap2 = [ voxmap( tmp2 <= cut-shift ), voxmap( tmp2 > cut-shift ) ];

lat_data = cnfield( mask, 1, voxmap2, 2, nsubj );

% approximate continuous field
params   = ConvFieldParams( FWHM * ones([1 D]), 1 );
params.lat_masked = lat_masked;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
d2cfield = convfield( lat_data, params, 2 );
LKC1     = LKC_voxmfd_est( cfield, dcfield, d2cfield );
LKC1_HPE = LKC_HP_est( cfield, Mboot, 1 );

params = ConvFieldParams( FWHM * ones([1 D]), 3 );
params.lat_masked = lat_masked;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
d2cfield = convfield( lat_data, params, 2 );
LKC3 = LKC_voxmfd_est( cfield, dcfield, d2cfield );
LKC3_HPE = LKC_HP_est( cfield, Mboot, 1 );

params = ConvFieldParams( FWHM * ones([1 D]), 5 );
params.lat_masked = lat_masked;
cfield   = convfield( lat_data, params, 0 );
dcfield  = convfield( lat_data, params, 1 );
d2cfield = convfield( lat_data, params, 2 );
LKC5 = LKC_voxmfd_est( cfield, dcfield, d2cfield );
LKC5_HPE = LKC_HP_est( cfield, Mboot, 1 );

% Here we see that the second integral somehow must be wrong
ConvE = [ theoryL; LKC1; LKC3; LKC5 ]'
HPE = [ theoryL; LKC1_HPE.hatL'; LKC3_HPE.hatL'; LKC5_HPE.hatL' ]'