%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the Reimmetric_est function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D stationary example
FWHM = 6; resadd = 1; nsubj = 100; nvox = 100;
params = ConvFieldParams(FWHM, resadd);
trunc = ceil(4*FWHM2sigma(FWHM));
noise = wfield(nvox + 2*trunc, nsubj);
smooth_f = convfield(noise, params);
start = floor((trunc/(nvox+2*trunc))*smooth_f.masksize(1));
stat_f = smooth_f(start:smooth_f.masksize(1)-start);
smooth_f_deriv = convfield(noise, params, 1);
stat_f_derivs = smooth_f_deriv(start:smooth_f.masksize(1) -start);

% Obtain the answer via Riemmetric_est
Lambda = mean(Riemmetric_est( stat_f, stat_f_derivs ).field);

% Obtain the answer from theory from the FWHM
FWHM2Lambda(FWHM,1)

%% %% 2D Examples
%% Simple 2D example


%% %% 3D Examples
%% Simple 3D example