%% Test the cope set code
%--------------------------------------------------------------------------
clear all
close all

%% Generate data
%--------------------------------------------------------------------------
% Dimension of the domain
D = 2;
% sample size
N = 200;
% Number of simulations
M = 100;
% size of domain
dim    = 100 * ones( [ 1 D ] );
% FWHM first smoothing
FWHM_1 = 3;
% FWHM second smoothing
FWHM_2 = 12;

% standard deviation of error
switch D
    case 1
        sd = 2;
        % voxelmap (identity)
        voxmap = 1:prod( dim );
    case 2
        sd = 10;
        % voxelmap (identity)
        voxmap = reshape( 1:prod( dim ), dim );
    case 3
        sd = 50;
        % voxelmap (identity)
        voxmap = reshape( 1:prod( dim ), dim );
end


% Sample grid
switch D
    case 1
        X = 1:dim(1);
    case 2
        [ X, Y ] = meshgrid( 1:dim(1), 1:dim(2) );
    case 3
        [ X, Y, Z ] = meshgrid( 1:dim(1), 1:dim(2), 1:dim(3) );
end

% population mean
mu = X / 25;
if D == 1
    mu = mu';
end

% create a smooth random field with linear mean
field = mu + sd * cnfield( dim, FWHM_1, voxmap, FWHM_2, N*M );
field = reshapeFiber(field, [N, M]);

%% Cope Set estimation
%--------------------------------------------------------------------------
% Threshold for Cope Set
c    = 1.5;
lvls = 0.90;

thresh    = zeros([field.masksize, 1, M, 2]);
quantiles = zeros([1 M]);
hatmu     = zeros([field.masksize, M]);
hatsigma  = zeros([field.masksize, M]);
len_bdry_tmp = zeros([1 M]);

% Get the Cope Set thresholds etc
for m = 1:M
    [ thresh_tmp, quantiles_tmp, hatmu_tmp, hatsigma_tmp, len_bdry_tmp ] = ...
                                                    CopeSets( field(:,:,:,m), c, lvls );
    thresh(:,:,1,m,:) = thresh_tmp;
    quantiles(m)    = quantiles_tmp;
    hatmu(:,:,m)    = hatmu_tmp.field;
    hatsigma(:,:,m) = hatsigma_tmp.field;
    len_bdry_tmp(m) = len_bdry_tmp;                                                
end

%% Quick and dirty plots
%--------------------------------------------------------------------------
switch D
    case 1
        % plot the mean
        figure(1), clf,
        subplot(2,2,1)
        plot( field )
        title( "Samples of observations" )
        subplot(2,2,2)
        plot( field- mu )
        title( "Samples of error field" )

        subplot(2,2,3)
        plot( hatmu )
        title( "estimated mean" )
        subplot(2,2,4)
        plot( hatsigma )
        title( "estimated std" )

        % Make a quick plot
        low_bd = hatmu.field > thresh(:,:,1);
        est    = ( hatmu.field > c );
        up_bd  = ( hatmu.field > thresh(:,:,2) );

        step_plot = low_bd + est + up_bd;

        figure(2)
        plot( step_plot )
        title( "Cope Sets: lower set = ( 1, 2, 3 ), est exc. set = ( 2, 3 ), upper set = ( 3 )." )
    case 2
        % plot the mean
        figure(1)
        subplot(2,2,1)
        imagesc( field( :, :, 1, 1 ) ), colorbar
        title( "sample of observations" )
        subplot(2,2,2)
        imagesc( field( :, :, 1, 1 ) - mu ), colorbar
        title( "Sample of error field" )

        subplot(2,2,3)
        imagesc( hatmu_tmp ), colorbar
        title( "estimated mean" )
        subplot(2,2,4)
        imagesc( hatsigma_tmp ), colorbar
        title( "estimated std" )

        % Make a quick plot
        low_bd = hatmu_tmp.field > thresh_tmp(:,:,1);
        est    = ( hatmu_tmp.field > c );
        up_bd  = ( hatmu_tmp.field > thresh_tmp(:,:,2) );

        step_plot = low_bd + est + up_bd;

        figure(2)
        imagesc( step_plot ), colorbar
        title( "Cope Sets: lower set = ( 1, 2, 3 ), est exc. set = ( 2, 3 ), upper set = ( 3 )." )
    case 3
        slice = 10;
        % plot the mean
        figure(1)
        subplot(2,2,1)
        imagesc( field( :, :, slice, 1 ) ), colorbar
        title( "sample of observations" )
        subplot(2,2,2)
        imagesc( field( :, :, slice, 1 ) - mu(:,:,slice) ), colorbar
        title( "Sample of error field" )

        subplot(2,2,3)
        imagesc( hatmu( :, :, slice ) ), colorbar
        title( "estimated mean" )
        subplot(2,2,4)
        imagesc( hatsigma( :, :, slice ) ), colorbar
        title( "estimated std" )

        % Make a quick plot
        low_bd = hatmu.field > thresh(:,:,1);
        est    = ( hatmu.field > c );
        up_bd  = ( hatmu.field > thresh(:,:,2) );

        step_plot = low_bd + est + up_bd;

        figure(2)
        imagesc( step_plot( :, :, slice ) ), colorbar
        title( "Cope Sets: lower set = ( 1, 2, 3 ), est exc. set = ( 2, 3 ), upper set = ( 3 )." )

end

%% Checking whether the CoPE set covers the true level set

% only voxels
coveringRateNaive = CovRateLvlSets( mu, hatmu, thresh, c, 0 )

% linear interpolated voxels
coveringRateInterpol = CovRateLvlSets( mu, hatmu, thresh, c,1 )

