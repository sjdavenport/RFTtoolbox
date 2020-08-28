%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the XXX function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D example
Dim = 10;
nsubj = 100;
noise = noisegen(Dim, nsubj, 6);
[ fwhm_est, Lambda_est, sigma_est, fwhm_est_forman ] = est_smooth(noise)

%%
Dim = [250,250];
nsubj = 1;
noise = noisegen(Dim, nsubj, 20);
[ fwhm_est, Lambda_est, sigma_est, fwhm_est_forman ] = est_smooth(noise)

%%
Dim = [250,250];
nsubj = 20;
noise = noisegen(Dim, nsubj, 4);
[ fwhm_est, Lambda_est, sigma_est, fwhm_est_forman ] = est_smooth(noise)

%%
Dim = [250,250];
nsubj = 100;
noise = noisegen(Dim, nsubj, 6);
est_smooth(noise)

Dim = [91,109,91];
nsubj = 100;
noise = noisegen(Dim, nsubj, 6);
est_smooth(noise) %Gets around 6.14

% 3D example with a mask
Dim = [91,109,91];
nsubj = 20;
noise = noisegen(Dim, nsubj, 3);
mask = imgload('MNImask');
est_smooth(noise, mask) %Gets around 3.26