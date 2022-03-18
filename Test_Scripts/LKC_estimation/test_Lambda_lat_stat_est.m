%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the Lambda_lat_stat_est function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 2D Examples

%% Simple 2D example with many subjects - comparing Kiebel estimators
dim = [100,100];
nsubj = 20;
lat_data = wfield(dim, nsubj);
FWHM = 3;
params = ConvFieldParams([FWHM, FWHM], 0);
f = convfield(lat_data, params);
[ Lambda ] = Lambda_lat_stat_est( f );
[ Lambda_onesided ] = Lambda_lat_stat_est( f, "onesided" );
FWHM_lat_stat_kiebel = sqrt(4*log(2)/mean(diag(Lambda)))
FWHM_lat_stat_onesided_kiebel = sqrt(4*log(2)/mean(diag(Lambda_onesided)))

[ ~, fwhm_est_smooth_kiebel] = est_smooth( f.field, 0);
fwhm_est_smooth_kiebel = mean(fwhm_est_smooth_kiebel)

% The one-sided estimates seems to be also identical. The two-sided
% estimates of Lambda_lat_stat_est seems to be quite off

%% Simple 2D example with 1 subject
dim = [100,100];
lat_data = wfield(dim);
FWHM = 3;
params = ConvFieldParams([FWHM, FWHM], 0);
f = convfield(lat_data, params);
[ Lambda ] = Lambda_lat_stat_est( f )
[ Lambda_onesided ] = Lambda_lat_stat_est( f, "onesided" )
FWHM_est = sqrt(4*log(2)/mean(diag(Lambda)))
FWHM_est_onesided = sqrt(4*log(2)/mean(diag(Lambda_onesided)))

[ fwhm_est_forman, fwhm_est_kiebel] = est_smooth( f.field, 0)

% I guess it's not set up for dealing with 1 subject yet!



%% %% 3D Examples
%% Simple 3D example