%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the ECcurve function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D example
nvox = 100; FWHM = 5; D = length(Dim); resadd = 1;
lat_data = wnfield(nvox);
params = ConvFieldParams(FWHM, resadd);
smooth_data = convfield(lat_data, params);
[ curve, x ] = ECcurve(smooth_data, [-5,5], 0.001);
[ curve_dep, x ] = ECcurve_dep(smooth_data, [-5,5], 0.001);

plot(x,curve_dep)
hold on
plot(x, curve)

%% %% 2D Examples
%% Simple 2D example
Dim = [10,10]; FWHM = 5; D = length(Dim);
lat_data = wnfield(Dim);
[ curve, x ] = ECcurve(lat_data);
plot(x, curve)


%% %% 3D Examples
%% Simple 3D example
Dim = [100,100,100]; FWHM = 5; D = length(Dim); resadd = 0;
lat_data = wnfield(Dim);
params = ConvFieldParams(repmat(FWHM,1,D), resadd);
smooth_data = convfield(lat_data, params);
limits = [-5,5]; increm = 0.001;
tic
[ curve, x ] = ECcurve(smooth_data, limits, increm);
toc
% plot(x,curve)