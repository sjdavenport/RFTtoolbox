%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests compares voxelwise permutation and RFT
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Simple 1D example (already with resadd = 1 RFT does better!)
signal = 2*[ zeros(1,10), ones(1,5), zeros(1,10) ]'; 
nvox = length(signal); nsubj = 75; FWHM = 2; 
noisey_data = 10*randn([nvox,nsubj]) + signal;

% RFT thresholding
[im, threshold, ~, xvals] = vRFT( noisey_data, FWHM );

% Permutation thresholidng
threshold_perm = perm_thresh(noisey_data, 'T', FWHM);
im_perm = mvtstat( noisey_data ) > threshold_perm;

% Plot results
subplot(3,1,1)
plot(xvals{1}, im); title('RFT Discoveries')
xlim(vec2lim(xvals{1}))
subplot(3,1,2)
plot(fconv(signal, FWHM)); title('Smoothed Signal');
xlim(vec2lim(xvals{1}))
subplot(3,1,3)
plot(im_perm); title('Permutation Discoveries')
