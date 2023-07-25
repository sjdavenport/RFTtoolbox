%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the tfce function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D example


%% %% 2D Examples
dim = [50,50]; nsubj = 50; FWHM = 0;
Sig = 0.25*peakgen(1, 10, 8, dim);
Sig = 0.5*square_signal(dim, 4, {[25,20], [25,30]} );
data = wfield(dim, nsubj);
data.field = data.field + Sig;
tstat = convfield_t(data, FWHM);
tstat_tfce = tfce(tstat.field,2,0.5,8,0.05);

subplot(1,2,1)
surf(tstat.field)
title('Original tstat')
view([-14,15])
subplot(1,2,2)
surf(tstat_tfce)
title('TFCE')
view([-14,15])
fullscreen
%% Smooth after
tstat_im = mvtstat(data.field);
smooth_tstat = fconv(tstat_im, FWHM, 2);
smooth_tstat_tfce = tfce(smooth_tstat,2,0.5,8,0.05);
subplot(1,2,1)
surf(smooth_tstat)
subplot(1,2,2)
surf(smooth_tstat_tfce)

%% %% 3D Examples
%% Simple 3D example