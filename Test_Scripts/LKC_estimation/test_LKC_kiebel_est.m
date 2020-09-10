%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the LKC_kiebel_est function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D examples
% Generate a random field
D = 1;
params   = ConvFieldParams( 6*ones( [ 1 D ] ), 1, 1, false );
lat_data = wfield( pad_vals( true( [ 25 1 ] ), 5, false ), 30  );
field    = convfield( lat_data, params );

% Estimate LKCs, Riemannian metric and FWHM using Kiebel 1999, Note that if
% we would subtract the mean the second input would be 1.
[ Lkiebel, Lambda, FWHM ] = LKC_kiebel_est( field, 0 );
LKC_forman_est( field, 0 )
% Estimate LKCs using convolution approach
Lconv = LKC_latconv_est( lat_data, params );

% Kiebel is slightly below convolution estimator all the time
[ Lkiebel, Lconv ] 

%% %% 2D example
% Generate a random field
D = 2;
params   = ConvFieldParams( 6 * ones( [ 1 D ] ), 1, 1, false );
lat_data = wfield( pad_vals( true( 25 * ones([ 1 D ]) ), 5, false ), 30  );
field    = convfield( lat_data, params );

% Estimate LKCs, Riemannian metric and FWHM using Kiebel 1999, Note that if
% we would subtract the mean the second input would be 1.
[ Lkiebel, Lambda, FWHM ] = LKC_kiebel_est( field, 0 );
[ Lforman, Lkiebel2 ] = LKC_forman_est( field, 0 );


% Estimate LKCs using convolution approach
Lconv = LKC_latconv_est( lat_data, params );

% mostly Kiebel is slightly below convolution estimator all the time
[ Lkiebel; Lconv ]
FWHM


%% %% 3D example
% Generate a random field
D = 3;
params   = ConvFieldParams( 6 * ones( [ 1 D ] ), 1, 1, false );
lat_data = wfield( pad_vals( true( 25 * ones([ 1 D ]) ), 5, false ), 30  );
field    = convfield( lat_data, params );

% Estimate LKCs, Riemannian metric and FWHM using Kiebel 1999, Note that if
% we would subtract the mean the second input would be 1.
[ Lkiebel, Lambda, FWHM ] = LKC_kiebel_est( field, 0 );
LKC_forman_est( field, 0 )

% Estimate LKCs using convolution approach
Lconv = LKC_latconv_est( lat_data, params );

% mostly Kiebel is slightly below convolution estimator all the time
[ Lkiebel; Lconv ]
