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
FWHM   = 3;
pad    = ceil( 4*FWHM2sigma( FWHM ) );

%% Example with recangular mask
% Generate mask
mask = pad_vals( true( [ T, 1 ] ), pad, false );
lat_masked = false;

% Generate params object for convfields
resThy = 301;
params = ConvFieldParams( FWHM, resThy, ceil(resThy/2), lat_masked );

% LKC from continuous theory
theoryLt = LKC_isogauss_theory( FWHM, T  ); % This is for continuous fields
theoryL  = LKC_wncfield_theory( mask, params );

% Generate test data
lat_data = wfield( mask, nsubj );

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
struct('theory', theoryL,...
       'theoryIso', theoryLt,...
       'res1', LKC1,...
       'res3', LKC3,...
       'res5', LKC5 )

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
lat_data = wfield( mask, nsubj );

% Estimate across different resadd
params = ConvFieldParams( FWHM, 0, 0, false );
LKC0   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( FWHM, 1, 0, false );
LKC1   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( FWHM, 2, 0, false );
LKC2   = LKC_latconv_est( lat_data, params );

% Values are stable accross different resadd increases.
% Note that resadd should be odd. Otherwise the underlying
% manifold changes
% Compare results
struct('theory', theoryL,...
       'theoryLt', theoryLt,...
       'res1', LKC0,...
       'res3', LKC1,...
       'res5', LKC2 )

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
lat_data = wfield( mask, nsubj );

% Estimate across different resadd
params = ConvFieldParams( FWHM, 1, ceil(1/2), false );
LKC1   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( FWHM, 3, ceil(3/2), false );
LKC3   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( FWHM, 5, ceil(5/2), false );
LKC5   = LKC_latconv_est( lat_data, params );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
struct('theory', theoryL,...
       'theoryIso', theoryLt,...
       'res1', LKC1,...
       'res3', LKC3,...
       'res5', LKC5 )

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
lat_data = wfield( mask, nsubj );

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
% should be odd.struct('theory', theoryL,...
struct('theory', theoryL,...
       'theoryIso', theoryLt,...
       'res1', LKC1,...
       'res3', LKC3,...
       'res5', LKC5 )

%% %% D = 2 
%% % Parameters for the field
T      = 49;
nsubj  = 100;
FWHM   = 3%sigma2FWHM(1.2);
pad = ceil( 4*FWHM2sigma( FWHM ) );

%% Example with recangular mask
% Get mask
T = 49
mask = pad_vals( true( [ T T ] ), pad, false );

% Mask the lattice before smoothing or not
mask_lat = false;
%mask_lat = true;

% Generate params object for convfields
params = ConvFieldParams( [ FWHM, FWHM ], 5, ceil(5/2), mask_lat );

% LKC from continuous theory
theoryLt = LKC_isogauss_theory( FWHM, [ T T ]  );
theoryL  = LKC_wncfield_theory( mask, params );

% Generate test data
lat_data = wfield( mask, nsubj );

% Estimate across different resadd
params = ConvFieldParams( [ FWHM, FWHM ], 1, ceil(1/2), mask_lat );
LKC1   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( [ FWHM, FWHM ], 3, ceil(3/2), mask_lat );
LKC3   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( [ FWHM, FWHM ], 5, ceil(5/2), mask_lat );
LKC5   = LKC_latconv_est( lat_data, params );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
struct('theoryIso', theoryLt,...
       'theoryWn', theoryL,...
        'res1', LKC1,...
        'res3', LKC3,...
        'res5', LKC5 )

%% Example with recangular mask no enlargement
% Get mask
mask = pad_vals( true( [ T T ] ), pad, false );

% Generate params object for convfields
params = ConvFieldParams( [ FWHM, FWHM ], 21, 0, false );

% LKC from continuous theory
theoryLt = LKC_isogauss_theory( FWHM, [ T T ]  );
theoryL  = LKC_wncfield_theory( mask, params );

% Generate test data
lat_data = wfield( mask, nsubj );

% Estimate across different resadd
params = ConvFieldParams( [ FWHM, FWHM ], 0, 0, false );
LKC0   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( [ FWHM, FWHM ], 1, 0, false );
LKC1   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( [ FWHM, FWHM ], 2, 0, false );
LKC2   = LKC_latconv_est( lat_data, params );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
struct('theoryIso', theoryLt,...
       'theoryWn', theoryL,...
        'res0', LKC0,...
        'res1', LKC1,...
        'res2', LKC2 )

%% Example with complicated mask
Sig =   peakgen([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
mask = logical( Sig > 0.02 & Sig < 1.1 );
figure(1), clf,
imagesc( mask ), colorbar,
title("mask")
clear Sig

lat_mask = true;

% Generate params object for convfields
params = ConvFieldParams( [ FWHM, FWHM ], 5, ceil(5/2), lat_mask );

% LKC from continuous theory
theoryL  = LKC_wncfield_theory( mask, params );

% Generate test data
lat_data = wfield( mask, nsubj );

% Estimate across different resadd
params = ConvFieldParams( [ FWHM, FWHM ], 1, ceil(1/2), lat_mask );
LKC1   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( [ FWHM, FWHM ], 3, ceil(3/2), lat_mask );
LKC3   = LKC_latconv_est( lat_data, params );
params = ConvFieldParams( [ FWHM, FWHM ], 5, ceil(5/2), lat_mask );
LKC5   = LKC_latconv_est( lat_data, params );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
struct('theoryIso', theoryLt,...
       'theoryWn', theoryL,...
        'res1', LKC1,...
        'res3', LKC3,...
        'res5', LKC5 )

    
%% %% D = 3 
% Parameters for the field
T      = 5;
nsubj  = 200;
FWHM   = sigma2FWHM(1.5);
FWHM   = 6;
pad    = ceil( 4*FWHM2sigma( FWHM ) );

nn = 40;

%% Rectangular domain example
mask = ones([ T T T ]);
mask = zeros([ T-1 T-1 T ]);

mask(2,2,2) = 1;
mask(2,2,3) = 1;
mask(2,2,4) = 1;
mask(2,3,4) = 1;
mask(3,3,2) = 1;
mask(3,3,4) = 1;

% Generate rectangular mask with a padded zero collar 
mask = pad_vals( mask , pad, false );
lat_masked = false;

% Get theoretical LKC
theory_res = 1;
params = ConvFieldParams( [ FWHM, FWHM, FWHM ], theory_res, ceil(theory_res/2), lat_masked );
[theoryL, g_thy]  = LKC_wncfield_theory( mask, params );
% LKC from continuous theory
theoryLt = LKC_isogauss_theory( FWHM, [ T T T ]  );

% Generate test data
lat_data = wfield( mask, nsubj );

% Estimate across different resadd
rr = 1;
params = ConvFieldParams( [ FWHM, FWHM, FWHM ], rr, ceil(rr/2), lat_masked );
[LKC1, ~, ~, ~,voxmfd_est]   = LKC_latconv_est( lat_data(:,:,:,1:nn), params );
params = ConvFieldParams( [ FWHM, FWHM, FWHM ], 3, ceil(3/2), lat_masked );
LKC3   = LKC_latconv_est( lat_data(:,:,:,1:nn), params );
params = ConvFieldParams( [ FWHM, FWHM, FWHM ], 5, ceil(5/2), lat_masked );
LKC5   = LKC_latconv_est( lat_data(:,:,:,1:nn), params );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
struct('theoryIso', theoryLt,...
       'theoryWn', theoryL,...
        'res1', LKC1,...
        'res3', LKC3,...
        'res5', LKC5 )

voxmfd_est = Mask(voxmfd_est);
squeeze(g_thy.field(26,25,25,:,:))
squeeze(voxmfd_est.g.field(26,25,25,:,:))

%% Sphere domain example
% Generate rectangular mask with a padded zero collar
mask = ones( [ T T T ] );
for t = 2:(T-1)
    mask(3:(T-1),2:(T-1),t) = 0;
end

nn = 20;

mask = pad_vals( mask, pad, false );
lat_masked = true;

% Get theoretical LKC
params = ConvFieldParams( [ FWHM, FWHM, FWHM ], 1, ceil(1/2), lat_masked );
theoryL  = LKC_wncfield_theory( mask, params );
% LKC from continuous theory
theoryLt = LKC_isogauss_theory( FWHM, [ T T T ]  );

% Generate test data
lat_data = wfield( mask, nsubj );
cfield  = convfield( lat_data, params );

% HPE L_1
HPE  = LKC_HP_est( cfield, 1, 1 );

% Estimate across different resadd
params = ConvFieldParams( [ FWHM, FWHM, FWHM ], 1, ceil(1/2), lat_masked );
LKC1   = LKC_latconv_est( lat_data(:,:,:,1:nn), params );
params = ConvFieldParams( [ FWHM, FWHM, FWHM ], 3, ceil(3/2), lat_masked );
LKC3   = LKC_latconv_est( lat_data(:,:,:,1:nn), params );
params = ConvFieldParams( [ FWHM, FWHM, FWHM ], 5, ceil(5/2), lat_masked );
LKC5   = LKC_latconv_est( lat_data(:,:,:,1:nn), params );

% Values are stable accross different resadd increases. Note that resadd
% should be odd.
struct('theoryIso', theoryLt,...
       'theoryWn', theoryL,...
       'HPE', HPE.hatL',...
        'res1', LKC1,...
        'res3', LKC3,...
        'res5', LKC5 )
    
    

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