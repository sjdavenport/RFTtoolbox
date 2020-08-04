%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the vRFT function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D examples
%% Simple 1D example
signal = 2*[ zeros(1,10), ones(1,5), zeros(1,10) ]';
nvox = length(signal); nsubj = 75; FWHM = 1; resadd = 3;
lat_data = 10*wnfield(nvox, nsubj) + signal;
params = ConvFieldParams( FWHM, resadd );

% Threshold
[im, threshold, maxvals, L] = vRFT( lat_data, params, 1 );

% Plot results
subplot(2,1,1)
plot(im.xvals{1}, im.field); title('Discoveries')
xlim(vec2lim(im.xvals{1}))
subplot(2,1,2)
plot(fconv(signal, FWHM)); title('Smoothed Signal');
xlim(vec2lim(im.xvals{1}))
threshold

%% 1D example with no peak finding
[im, threshold, maxvals, L] = vRFT( lat_data, params, 0);

%% 1D example with mask
signal = 2*[ zeros(1,10), ones(1,5), zeros(1,10), ones(1,5), zeros(1,10)]'; 
nvox = length(signal); nsubj = 75; FWHM = 3; 


%% %% 2D examples
% Signal generation
signal_magnitude = 2; signal_radii = 3; smoothnessofpeaks = 10; 
Dim = [50,50]; peaklocation = {[25,25]};
signal = gensig( signal_magnitude, signal_radii, smoothnessofpeaks, ...
    Dim, peaklocation);

% Data generation and thresholding
nsubj = 75; resadd = 1;
noisey_data = 10*wnfield(Dim, nsubj) + signal;
FWHM = 3; smoothed_data = fconv(noisey_data.field, FWHM);

params = ConvFieldParams( [FWHM, FWHM], resadd );
[im, threshold, ~, L] = vRFT( noisey_data, params );

% Plot data and thresholded image
clf; subplot(4,1,1)
surf(signal); title('True Signal')
subplot(4,1,2)
surf(fconv(noisey_data.field(:,:,1),FWHM)); title('Single Subject')
subplot(4,1,3)
surf(convfield_t(noisey_data, params).field*(1/sqrt(nsubj)))
title('one sample Cohen''s d = T/nsubj^{1/2}')
subplot(4,1,4)
surf(im.field); title('Voxelwise RFT Activation')