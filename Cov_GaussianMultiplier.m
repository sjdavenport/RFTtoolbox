%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%    Testing the multiplier bootstrap process
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
% Get the toolboxes
addpath(genpath("/home/drtea/matlabToolboxes/BrainStat/"))
addpath(genpath("/home/drtea/matlabToolboxes/RFTtoolbox/"))
addpath(genpath("/home/drtea/matlabToolboxes/spm12/"))

% Load the data
data_path = '/home/drtea/Desktop/';
load([data_path 'UKB_2D_randomized'])

% Load MNImask
MNImask = imgload('/home/drtea/matlabToolboxes/BrainStat/BrainImages/MNImask');
MNImask_2D = MNImask(:,:,45);
mask = dilate_mask(MNImask_2D, -2);
figure
subplot(1,2,1)
imagesc(mask)
subplot(1,2,2)
imagesc(MNImask_2D)

if exist('im_store_2D', 'var')
    im_store = im_store_2D;
end

%% %% 1D
% Parameters simulation
D = 1;
niters = 500;
FWHM   = 3;
sample_size = 100;
resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );

% Get sliced data
data1D        = squeeze( im_store( :, 30, : ) );
mask1D        = squeeze( MNImask_2D( :, 30 ) );
mask1D_eroded = squeeze( mask( :, 30 ) );

% Summary stats
data_m = mean( data1D, 2 );
data_s = std( data1D, 0, 2 );

% normalized data
data1D_normalized = (data1D - data_m ) ./ data_s;
figure
plot( data1D )
figure
plot( mask1D )

%% Get base fields for the simulations
spfn_orig = get_sample_fields( data1D_normalized, mask1D_eroded, D );
base_fields = spfn_orig(20);

%% %% LKC estimation from the data
% Get a sample of the multiplier bootstrap field
spfn = @(nsubj)  Mask( multiplier_field( base_fields.lat_data, nsubj ) );
mfields = spfn(100);

[ L, L0, ~, cfields ] = LKC_latconv_est( mfields, params );

L_HPE = LKC_HP_est( cfields{1}, 1000, 1 );

% (This seems to work)
[L, L_HPE.hatL ]

%% Gaussian analysis (Multiplier field)
coverage1D = record_coverage( spfn, sample_size, params, niters, 10)

% nsubj=50, niters=5000
%    conv: 0.0428
%     lat: 0.0382
% finelat: 0.0418

% nsubj=200, niters=5000
%    conv: 0.0394
%     lat: 0.0348
% finelat: 0.0382

% randomized
% nsubj=100, niters=50
%    conv: 0.0360
%     lat: 0.0320
% finelat: 0.0360


%% EC curve analysis
LKC_estimate = mean( coverage1D.storeLKCs, 2 )';
[ curve, x ] = maxECcurve( coverage1D.convmaxima, 0.1 );
figure
plot(x, curve);
curve_conv = EEC( x, LKC_estimate, L0, 'T', sample_size-1 );
hold on 
plot( x, curve_conv );
hold off
title('Multiplier field: EEC and maxima distribution')

%% Plot the number of maxima above the threshold (not that maximal 3 are
% reported)
figure
plot(coverage1D.maxabovethreshold)
title('Multiplier field: number of peaks above the threshold')

sum(coverage1D.maxabovethreshold) / niters

%% %% Gaussian analysis (Multiplier field all basis functions)
base_fields = Field(mask1D_eroded);
base_fields.field = data1D_normalized;
% Get a sample of the multiplier bootstrap field
spfn = @(nsubj)  Mask( multiplier_field( base_fields, nsubj ) );

coverage1Dall = record_coverage( spfn, sample_size, params, niters, 10)

% Seems to work way better than small samples, i.e. covariance structure
% seems to admit reasonable approximation

% nsubj=50, niters=5000
%    conv: 0.0508
%     lat: 0.0460
% finelat: 0.0496
% nsubj=200, niters=5000
%    conv: 0.0468
%     lat: 0.0422
% finelat: 0.0462

% randomized
% nsubj=100, niters=50
%    conv: 0.0480
%     lat: 0.0460
% finelat: 0.0480

%% EC curve analysis
LKC_estimate = mean( coverage1Dall.storeLKCs, 2 )';
[ curve, x ] = maxECcurve( coverage1Dall.convmaxima, 0.1 );
figure
plot(x, curve);
curve_conv = EEC( x, LKC_estimate, L0, 'T', sample_size-1 );
hold on 
plot( x, curve_conv );
hold off
title('Multiplier field: EEC and maxima distribution')

%% Plot the number of maxima above the threshold (not that maximal 3 are
% reported)
figure
plot(coverage1Dall.maxabovethreshold)
title('Multiplier field: number of peaks above the threshold')

sum(coverage1Dall.maxabovethreshold) / niters

%% White noise
spfn = @(nsubj)  Mask( wnfield( mask1D_eroded, nsubj ) );
coverage1Dwn = record_coverage( spfn, sample_size, params, niters)

% nsubj=50, niters=5000
%    conv: 0.0528
%     lat: 0.0464
% finelat: 0.0514
% nsubj=200, niters=5000
%    conv: 0.0512
%     lat: 0.0452
% finelat: 0.0496

%% EC curve analysis
LKC_estimate = mean( coverage1Dwn.storeLKCs, 2 )';
[ curve, x ] = maxECcurve( coverage1Dwn.convmaxima, 0.1 );
figure
plot(x, curve);
curve_conv = EEC( x, LKC_estimate, L0, 'T', sample_size-1 );
hold on 
plot( x, curve_conv );
hold off
title('White noise: EEC and maxima distribution')

%% Plot the number of maxima above the threshold (not that maximal 3 are
% reported)
figure
plot(coverage1Dwn.maxabovethreshold)
title('White Noise: number of peaks above the threshold')

sum(coverage1Dwn.maxabovethreshold)/niters

%% %% 2D
% Parameters simulation
D = 2;
niters = 500;
FWHM   = 3;
sample_size = 100;
resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );

% Get sliced data
data2D        = im_store;
mask2D        = MNImask_2D;
mask2D_eroded = mask;

% Summary stats
data_m = mean( data2D, 3 );
data_s = std( data2D, 0, 3 );

data2D_normalized = ( data2D - data_m ) ./ data_s;

figure
imagesc( data2D(:,:,1) )
figure
imagesc( mask2D )

spfn_orig   = get_sample_fields( data2D_normalized, mask2D_eroded, D );
base_fields = spfn_orig(100);

data2df = base_fields.lat_data;
data2df.mask = mask2D_eroded;

[ L, L0, ~, cfields ] = LKC_latconv_est( data2df, params );

L_HPE = LKC_HP_est( cfields{1}, 1000, 1 );

% (This seems to work)
[ L; L_HPE.hatL' ]

%% Get base fields for the simulations
spfn_orig   = get_sample_fields( data2D_normalized, mask2D_eroded, D );
base_fields = spfn_orig(20);

%% %% LKC estimation from the data
% Get a sample of the multiplier bootstrap field
spfn = @(nsubj)  Mask( multiplier_field( base_fields.lat_data, nsubj ) );
mfields = spfn(100);

[ L, L0, ~, cfields ] = LKC_latconv_est( mfields, params );

L_HPE = LKC_HP_est( cfields{1}, 1000, 1 );

% (This seems to work)
[ L; L_HPE.hatL' ]

%% Gaussian analysis (Multiplier field)
coverage2D = record_coverage( spfn, sample_size, params, niters, 10 )

% nsubj=50, niters=5000
%    conv: 0.0270
%     lat: 0.0192
% finelat: 0.0254
% nsubj=200, niters=5000
%    conv: 0.0282
%     lat: 0.0216
% finelat: 0.0264

%% EC curve analysis
LKC_estimate = mean( coverage2D.storeLKCs, 2 )';
[ curve, x ] = maxECcurve( coverage2D.convmaxima, 0.1 );
figure
plot(x, curve);
curve_conv = EEC( x, LKC_estimate, L0, 'T', sample_size-1 );
hold on 
plot( x, curve_conv );
hold off
title('Multiplier field: EEC and maxima distribution')

%% Plot the number of maxima above the threshold (not that maximal 3 are
% reported)
% This is part of the problem. There are quite frequently more than one
% local maxima above the threshold. This means that EEC is larger at that
% threshold than distribution of the maximum, which is what we observe.
% Compare to the white noise simulation later.
figure
plot(coverage2D.maxabovethreshold)
title('Multiplier field: number of peaks above the threshold')

sum(coverage2D.maxabovethreshold)/niters

%% %% Gaussian analysis (Multiplier field all basis functions)
base_fields = Field(mask2D_eroded);
base_fields.field = data2D_normalized;
% Get a sample of the multiplier bootstrap field
spfn = @(nsubj)  Mask( multiplier_field( base_fields, nsubj ) );

coverage2Dall = record_coverage( spfn, sample_size, params, niters, 10)

% Seems to work way better than small samples, i.e. covariance structure
% seems to admit reasonable approximation

% nsubj=50, niters=5000
%    conv: 0.0492
%     lat: 0.0378
% finelat: 0.0456
% nsubj=200, niters=5000
%    conv: 0.0446
%     lat: 0.0364
% finelat: 0.0426

% randomized
% nsubj=100, niters=500
%    conv: 0.0380
%     lat: 0.0280
% finelat: 0.0360

%% EC curve analysis
LKC_estimate = mean( coverage2Dall.storeLKCs, 2 )';
[ curve, x ] = maxECcurve( coverage2Dall.convmaxima, 0.1 );
figure
plot(x, curve);
curve_conv = EEC( x, LKC_estimate, L0, 'T', sample_size-1 );
hold on 
plot( x, curve_conv );
hold off
title('Multiplier field: EEC and maxima distribution')

%% Plot the number of maxima above the threshold (not that maximal 3 are
% reported)
figure
plot(coverage2Dall.maxabovethreshold)
title('Multiplier field: number of peaks above the threshold')

sum(coverage2Dall.maxabovethreshold) / niters


%% %% White noise
spfn = @(nsubj)  Mask( wnfield( mask2D_eroded, nsubj ) );
coverage2Dwn = record_coverage( spfn, sample_size, params, niters)

% nsubj=50, niters=5000
%    conv: 0.0494
%     lat: 0.0282
% finelat: 0.0444
% nsubj=200, niters=5000
%    conv: 0.0520
%     lat: 0.0308
% finelat: 0.0454

%% EC curve analysis
LKC_estimate = mean( coverage2Dwn.storeLKCs, 2 )';
[ curve, x ] = maxECcurve( coverage2Dwn.convmaxima, 0.1 );
figure
plot(x, curve);
curve_conv = EEC( x, LKC_estimate, L0, 'T', sample_size-1 );
hold on 
plot( x, curve_conv );
title('White noise: EEC and maxima distribution')
hold off

%% Plot the number of maxima above the threshold (not that maximal 3 are
% reported)
figure
plot( coverage2Dwn.maxabovethreshold )
title('White noise: number of peaks above the threshold')

sum(coverage2Dwn.maxabovethreshold)/niters