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
%% Simple 1D example single peak
nvox = 10; nsubj = 100; lat_data = wfield([nvox,1], nsubj);
peakspec = {[3,7]}; peakparams = {[2,2]};
[~, meanfn] = peakgen1D( 1:0.1:nvox, peakspec, peakparams, 1, 0.01);
FWHM = 2; out = convCR(lat_data, FWHM, {5}, meanfn);

%% Simple 1D example single peak multiple FWHM
nvox = 10; nsubj = 100; lat_data = wfield([nvox,1], nsubj);
peakspec = {[3,7]}; peakparams = {[2,2]};
[~, meanfn] = peakgen1D( 1:0.1:nvox, peakspec, peakparams, 1, 0.01);
FWHM = [2,3]; out = convCR(lat_data, [2,3], {5}, meanfn);

% %% Simple 1D example two peaks
% nvox = 15; nsubj = 100; lat_data = wfield([nvox,1], nsubj);
% peakspec = {[3,4], [4,5]}; peakparams = {[3,3]};
% [sig_store, meanfn] = peakgen1D( 1:0.1:nvox, peakspec, peakparams, 1, 0.01);
% plot(sig_store)
% FWHM = 2; 
% load([PIloc,'Variance/storevars'], 'allvars')
% lat_data.field = lat_data.field/sqrt(allvars(FWHM));
% out = convCR(lat_data, FWHM, {mean(peakspec{1}), mean(peakspec{2})}, meanfn);

%% %% 2D Examples
%% 2D example - one peak
nsubj = 50; dim = [30,30];
lat_data = wfield(dim, nsubj);
unsmooth_sig = peakgen(0.25, 1, 4, dim);
lat_data.field = lat_data.field + unsmooth_sig;
FWHM = 2; out2 = convCR(lat_data, FWHM, {[15,15]'})
subplot(2,2,1); surf(mean(lat_data).field)
subplot(2,2,2); surf(mean(convfield(lat_data, 2)).field)
subplot(2,2,3); surf(mean(convfield(lat_data, 4)).field)
subplot(2,2,4); surf(mean(convfield(lat_data, 6)).field)

%% Simple 2D example - two peaks
nsubj = 50; dim = [30,60];
lat_data = wfield(dim, nsubj);
unsmooth_sig = peakgen(1, 1, 4, dim, {[15,15],[15,45]});
lat_data.field = lat_data.field + unsmooth_sig;
FWHM = 2; out = convCR(lat_data, FWHM, {[15,15]', [15,45]'})
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
