%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the convpeakcov function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D example single peak
nvox = 10; nsubj = 100; lat_data = wfield([nvox,1], nsubj);
peakspec = {[3,7]}; peakparams = {[2,2]};
[~, meanfn] = peakgen1D( 1:0.1:nvox, peakspec, peakparams, 1, 0.01);
FWHM_vec = [2,3]; out = convpeakcov(lat_data, FWHM_vec, {5}, meanfn);

%%
out2 = convpeakcov(lat_data, 2, {5}, meanfn);
out4 = convpeakcov(lat_data, 4, {5}, meanfn);

%% %% 2D Examples
%% 2D example - one peak
nsubj = 50; dim = [30,30];
lat_data = wfield(dim, nsubj);
unsmooth_sig = peakgen(0.25, 1, 4, dim);
lat_data.field = lat_data.field + unsmooth_sig;
FWHM_vec = [2,3]; out = convpeakcov(lat_data, FWHM_vec, {[15,15]'});

%% Simple 2D example - two peaks
nsubj = 50; dim = [30,60];
lat_data = wfield(dim, nsubj);
unsmooth_sig = peakgen(1, 1, 4, dim, {[15,15],[15,45]});
lat_data.field = lat_data.field + unsmooth_sig;
FWHM = 2; out = convCR(lat_data, [2,3], {[15,15]', [15,45]'})
subplot(2,2,1); surf(mean(lat_data).field)
subplot(2,2,2); surf(mean(convfield(lat_data, 2)).field)
subplot(2,2,3); surf(mean(convfield(lat_data, 4)).field)
subplot(2,2,4); surf(mean(convfield(lat_data, 6)).field)

%% %% 3D Examples
%% Simple 3D example
nsubj = 50; dim = [30,30,30];
lat_data = wfield(dim, nsubj);
unsmooth_sig = peakgen(0.25, 1, 4, dim);
lat_data.field = lat_data.field + unsmooth_sig;
FWHM = 2; out = convCR(lat_data, FWHM, {[15,15,15]'})
subplot(2,2,1); surf(mean(lat_data).field)
subplot(2,2,2); surf(mean(convfield(lat_data, 2)).field)
subplot(2,2,3); surf(mean(convfield(lat_data, 4)).field)
subplot(2,2,4); surf(mean(convfield(lat_data, 6)).field)
