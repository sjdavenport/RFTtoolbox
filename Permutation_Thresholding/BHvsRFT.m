%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests compares power between voxelwise RFT and FDR
%%%    control using the BH procedure. Note that this comparison is not
%%%    entirely fair as FDR provides weaker control than FWER. However it's
%%%    still of interest. RFT should be more powerful when there isn't much
%%%    signal, despite being able to provie more stringent control!
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simple 1D example
% Note that RFT seems to be more powerful as the smoothness increases!
signal = 2*[ zeros(1,10), ones(1,5), zeros(1,10) ]'; mask = ones(length(signal),1);
nvox = length(signal); nsubj = 75; FWHM = 2; resadd = 5;
noisey_data = 10*randn([nvox,nsubj]) + signal;

% RFT thresholding
[im, threshold, ~, xvals] = vRFT( noisey_data, FWHM, mask, resadd );

% FDR thresholding
rej_locs = spatialBH( noisey_data, FWHM );

% Plot results
subplot(3,1,1)
plot(xvals{1}, im); title('RFT Discoveries')
xlim(vec2lim(xvals{1}))
subplot(3,1,2)
plot(fconv(signal, FWHM)); title('Smoothed Signal');
xlim(vec2lim(xvals{1}))
subplot(3,1,3)
plot(rej_locs); title('BH Discoveries')

%% Simple 1D example (already with resadd = 1 RFT does better!)
signal = 2*[ zeros(1,10), ones(1,5), zeros(1,10), ones(1,5) ]'; mask = ones(length(signal),1);
nvox = length(signal); nsubj = 75; FWHM = 2;
noisey_data = 10*randn([nvox,nsubj]) + signal;

% RFT thresholding
[im, threshold, ~, xvals] = vRFT( noisey_data, FWHM, mask, 3 );

% FDR thresholding
rej_locs = spatialBH( noisey_data, FWHM );

% Plot results
subplot(3,1,1)
plot(xvals{1}, im); title('RFT Discoveries')
xlim(vec2lim(xvals{1}))
subplot(3,1,2)
plot(fconv(signal, FWHM)); title('Smoothed Signal');
xlim(vec2lim(xvals{1}))
subplot(3,1,3)
plot(rej_locs); title('BH Discoveries')