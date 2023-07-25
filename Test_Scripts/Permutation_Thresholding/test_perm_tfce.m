%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the perm_tfce function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D example


%% %% 2D Examples
%% Simple 2D example
dim = [50,50]; nsubj = 50; FWHM = 0;
Sig = 0.25*peakgen(1, 10, 8, dim);
data = wfield(dim, nsubj).field + Sig;
threshold_tfce = perm_tfce(data, ones(dim));
tfce_tstat = tfce(data);
threshold_orig = perm_tfce(data, ones(dim));

%% %% 3D Examples
%% Simple 3D example