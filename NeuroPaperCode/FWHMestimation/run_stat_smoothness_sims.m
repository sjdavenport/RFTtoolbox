%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script runs the stationary smoothness simulations
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D 
nvox = 10; mask = true([nvox, 1]); nsubj = 100; FWHM_vec = 0.5:0.5:6;
global RFTboxloc
savedir = [RFTboxloc, 'NeuroPaperCode/FWHMestimation/oneDstatsims_nsubj', num2str(nsubj)];
stat_smoothness_sims( mask, nsubj, FWHM_vec, 1000, savedir )

%% %% 2D 
Dim = [10,10]; mask = true(Dim); nsubj = 10; FWHM_vec = 0.5:0.5:6;
savedir = [RFTboxloc, 'NeuroPaperCode/FWHMestimation/twoDstatsims'];
stat_smoothness_sims( mask, nsubj, FWHM_vec, 1000, savedir )

%% %% 3D 
Dim = [10,10,10]; mask = true([nvox, 1]); nsubj = 10; FWHM_vec = 0.5:0.5:6;
savedir = [RFTboxloc, 'NeuroPaperCode/FWHMestimation/threeDstatsims'];
stat_smoothness_sims( mask, nsubj, FWHM_vec, 1000, savedir )