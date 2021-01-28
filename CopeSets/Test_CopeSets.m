%% Test the cope set code
%--------------------------------------------------------------------------
clear all
close all

%% Generate data
%--------------------------------------------------------------------------
% sample size
N = 10;
% dimension of domain
dim    = [ 100 100 ];
% FWHM first smoothing
FWHM_1 = 3;
% FWHM second smoothing
FWHM_2 = 12;
% voxelmap (identity)
voxmap = reshape( 1:prod( dim ), dim );

% standard deviation of error
sd = 10;

% Sample grid
[ X, Y ] = meshgrid( 1:dim(1), 1:dim(2) );

% population mean
mu = X / 25;

% create a smooth random field with linear mean
field = mu + sd * cnfield( dim, FWHM_1, voxmap, FWHM_2, N );


%% Cope Set estimation
%--------------------------------------------------------------------------
% Threshold for Cope Set
c    = 1.2;
lvls = 0.90;

% Get the Cope Set thresholds etc
[ thresh, quantiles, hatdelta, hatsigma, len_bdry ] = CopeSets( field, c, lvls );

% plot the mean
figure(1)
subplot(2,2,1)
imagesc( field( :, :, 1 ) ), colorbar
title( "sample of observations" )
subplot(2,2,2)
imagesc( field( :, :, 1 ) - mu ), colorbar
title( "Sample of error field" )

subplot(2,2,3)
imagesc( hatdelta ), colorbar
title( "estimated mean" )
subplot(2,2,4)
imagesc( hatsigma ), colorbar
title( "estimated std" )


%% Quick and dirty plot
%--------------------------------------------------------------------------
% Make a quick plot
low_bd = hatdelta.field > thresh(:,:,1);
est    = ( hatdelta.field > c );
up_bd  = ( hatdelta.field > thresh(:,:,2) );

step_plot = low_bd + est + up_bd;

figure(2)
imagesc( step_plot ), colorbar
title( "Cope Sets: lower set = ( 1, 2, 3 ), est exc. set = ( 2, 3 ), upper set = ( 3 )." )