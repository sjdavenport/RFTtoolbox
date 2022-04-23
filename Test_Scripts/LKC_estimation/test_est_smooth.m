%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the est_smooth function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D example
Dim = 100;
nsubj = 100;
noise = noisegen(Dim, nsubj, 6);
[ fwhm_est, fwhm_est_forman ] = est_smooth(noise)

%% Add smoothing to detect smoothing
Dim = 100;
nsubj = 100;
FWHM = 4;
noise = noisegen(Dim, nsubj, FWHM);
[ fwhm_est, fwhm_est_forman ] = est_smooth(noise)

%% Single subject
Dim = [250,250];
nsubj = 1;
noise = noisegen(Dim, nsubj, 20);
[ fwhm_est, fwhm_est_forman ] = est_smooth(noise)

%% %% 2D Examples
%% Simple 2D example
Dim = [50,50];
nsubj = 10;
noise = noisegen(Dim, nsubj, 12);
[ fwhm_est, fwhm_est_forman ] = est_smooth(noise)

%% Add smoothing to detect smoothing
Dim = [50,50];
nsubj = 20;
noise = noisegen(Dim, nsubj, 4);
[ fwhm_est, fwhm_est_forman ] = est_smooth(noise)

%%
Dim = [250,250];
nsubj = 100;
noise = noisegen(Dim, nsubj, 6);
est_smooth(noise)

Dim = [91,109,91];
nsubj = 100;
noise = noisegen(Dim, nsubj, 6);
est_smooth(noise) %Gets around 6.14

%% 3D example with a mask
Dim = [91,109,91];
nsubj = 20;
noise = noisegen(Dim, nsubj, 3);
mask = imgload('MNImask');
est_smooth(noise, 0, mask) %Gets around 3.26