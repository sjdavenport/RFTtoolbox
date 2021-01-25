%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the stat_Lambda function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D example
FWHM = 10; nvox = 100; resadd = 0; nsubj = 100;
params = ConvFieldParams(FWHM, resadd);
lat_data = wfield([nvox, 1], nsubj);
cfields = convfield(lat_data, params);
dfields = convfield(lat_data, params, 1);
[Lambda, FWHM_est] = stat_Lambda( cfields, dfields )

%% %% 2D Examples
%% Simple 2D example
FWHM = 3; Dim = [100,100]; resadd = 0; nsubj = 100;
params = ConvFieldParams([FWHM,FWHM], resadd);
lat_data = wfield(Dim, nsubj);
cfields = convfield(lat_data, params);
dfields = convfield(lat_data, params, 1);
[Lambda, FWHM_est] = stat_Lambda( cfields, dfields );

det(Lambda)
(4*log(2))^2/FWHM^4

% Note that the edge effect means that the estimate of the FWHM is biased
% high!
%% %% 3D Examples
%% Simple 3D example
FWHM = 3; Dim = [50,50,50]; resadd = 0; nsubj = 100;
params = ConvFieldParams([FWHM,FWHM, FWHM], resadd);
lat_data = wfield(Dim, nsubj);
cfields = convfield(lat_data, params);
dfields = convfield(lat_data, params, 1);
[Lambda, FWHM_est] = stat_Lambda( cfields, dfields );

det(Lambda)
(4*log(2))^3/FWHM^6