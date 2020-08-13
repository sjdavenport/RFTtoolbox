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
load([data_path 'UKB_2D'])

% Load MNImask
MNImask = imgload('/home/drtea/matlabToolboxes/BrainStat/BrainImages/MNImask');
MNImask_2D = MNImask(:,:,45);
mask = dilate_mask(MNImask_2D, -2);
figure
subplot(1,2,1)
imagesc(mask)
subplot(1,2,2)
imagesc(MNImask_2D)

%% %% 1D
% Parameters simulation
D = 1;
niters = 500;
FWHM   = 3;
sample_size = 30;
resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );

% Get sliced data
data1D        = squeeze( im_store( :, 30, : ) );
mask1D        = squeeze( MNImask_2D( :, 30 ) );
mask1D_eroded = squeeze( mask( :, 30 ) );

% Summary stats
data_m = mean( data1D, 2 );
data_s = std( data1D, 0, 2 );

figure
plot( data1D )
figure
plot( mask1D )

%% Get base fields for the simulations
spfn_orig = get_sample_fields( data1D, mask1D_eroded, D );
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
coverage1D = record_coverage( spfn, sample_size, params, niters)

% nsubj=300, niters=1000
% coverage1D = 
% 
%   struct with fields:
% 
%              conv: 0.0390
%               lat: 0.0370
%           finelat: 0.0390
%     finelatmaxima: [1×1000 double]
%         latmaxima: [1×1000 double]
%        convmaxima: [1×1000 double]
%         storeLKCs: [1×1000 double]

%% EC curve analysis
LKC_estimate = mean( coverage1D.storeLKCs, 2 )';
[ curve, x ] = maxECcurve( coverage1D.convmaxima, 0.1 );
plot(x, curve);
curve_conv = EEC( x, LKC_estimate, L0, 'T', sample_size-1 );
hold on 
plot( x, curve_conv );
hold off
title('Multiplier field: EEC and maxima distribution')

%% Plot the number of maxima above the threshold (not that maximal 3 are
% reported)
tmp = coverage1D.allmaxima >= repmat( coverage1D.thresholds, [ 3 1 ] );
figure
plot(sum(tmp))
title('Multiplier field: number of peaks above the threshold')

tmp = coverage1D.convmaxima >= coverage1D.thresholds;
sum(tmp)/niters

%% White noise
spfn = @(nsubj)  Mask( wnfield( mask1D_eroded, nsubj ) );
coverage1D = record_coverage( spfn, sample_size, params, niters)

% nsubj=300, niters=1000
% coverage1D = 
% 
%   struct with fields:
% 
%              conv: 0.0390
%               lat: 0.0370
%           finelat: 0.0390
%     finelatmaxima: [1×1000 double]
%         latmaxima: [1×1000 double]
%        convmaxima: [1×1000 double]
%         storeLKCs: [1×1000 double]

%% EC curve analysis
LKC_estimate = mean( coverage1D.storeLKCs, 2 )';
[ curve, x ] = maxECcurve( coverage1D.convmaxima, 0.1 );
plot(x, curve);
curve_conv = EEC( x, LKC_estimate, L0, 'T', sample_size-1 );
hold on 
plot( x, curve_conv );
hold off
title('White noise: EEC and maxima distribution')

%% Plot the number of maxima above the threshold (not that maximal 3 are
% reported)
tmp = coverage1D.allmaxima >= repmat( coverage1D.thresholds, [ 3 1 ] );
figure
plot(sum(tmp))
title('White Noise: number of peaks above the threshold')


tmp = coverage1D.convmaxima >= coverage1D.thresholds;
sum(tmp)/niters


%% %% 2D
close all
% Parameters simulation
D = 2;
niters = 500;
FWHM   = 3;
sample_size = 20;
resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );

% Get sliced data
data2D        = im_store;
mask2D        = MNImask_2D;
mask2D_eroded = mask;

% Summary stats
data_m = mean( data2D, 3 );
data_s = std( data2D, 0, 3 );

figure
imagesc( data2D(:,:,1) )
figure
imagesc( mask2D )

%% Get base fields for the simulations
spfn_orig   = get_sample_fields( data2D, mask2D_eroded, D );
base_fields = spfn_orig(20);

%%
data_path = 'C:\Users\12SDa\davenpor\Data\RestingStateData\EklundData\';
load([data_path 'UKB_2D_randomized'])
MNImask = imgload('MNImask');
MNImask_2D = MNImask(:,:,45);

D = 2; niters = 1000;
spfn = get_sample_fields( im_store_2D, MNImask_2D, D );
base_fields = spfn(20);

%% %% LKC estimation from the data
% Get a sample of the multiplier bootstrap field
spfn = @(nsubj)  Mask( multiplier_field( base_fields.lat_data, nsubj ) );
mfields = spfn(100);

[ L, L0, ~, cfields ] = LKC_latconv_est( mfields, params );

L_HPE = LKC_HP_est( cfields{1}, 1000, 1 );

% (This seems to work)
[ L; L_HPE.hatL' ]

%% Gaussian analysis (Multiplier field)
coverage2D = record_coverage( spfn, sample_size, params, niters, 10)

% sample_size = 20, niters = 1000
% coverage2D = 
%   struct with fields:
% 
%              conv: 0.0390
%               lat: 0.0230
%           finelat: 0.0330
%     finelatmaxima: [1×1000 double]
%         latmaxima: [1×1000 double]
%        convmaxima: [1×1000 double]
%         storeLKCs: [2×1000 double]

%%
ECcurveanal( coverage2D, MNImask_2D, sample_size)

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
tmp = coverage2D.allmaxima >= repmat( coverage2D.thresholds, [ 3 1 ] );

% This is part of the problem. There are quite frequently more than one
% local maxima above the threshold. This means that EEC is larger at that
% threshold than distribution of the maximum, which is what we observe.
% Compare to the white noise simulation later.
figure
plot(sum(tmp))
title('Multiplier field: number of peaks above the threshold')

tmp = coverage2D.convmaxima >= coverage2D.thresholds;
sum(tmp)/niters

%% %% White noise
spfn = @(nsubj)  Mask( wnfield( mask2D_eroded, nsubj ) );
coverage2D = record_coverage( spfn, sample_size, params, niters)

%% EC curve analysis
LKC_estimate = mean( coverage2D.storeLKCs, 2 )';
[ curve, x ] = maxECcurve( coverage2D.convmaxima, 0.1 );
figure
plot(x, curve);
curve_conv = EEC( x, LKC_estimate, L0, 'T', sample_size-1 );
hold on 
plot( x, curve_conv );
title('White noise: EEC and maxima distribution')
hold off

%% Plot the number of maxima above the threshold (not that maximal 3 are
% reported)
tmp = coverage2D.allmaxima >= repmat( coverage2D.thresholds, [ 3 1 ] );
figure
plot(sum(tmp))
title('White noise: number of peaks above the threshold')


tmp = coverage2D.convmaxima >= coverage2D.thresholds;
sum(tmp)/niters