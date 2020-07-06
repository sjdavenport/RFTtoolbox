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
threshold_perm = perm_thresh(noisey_data, 'T', FWHM, NaN, 1);
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

%% Test Power (1D)
niters = 1000;
signal = 2*[ zeros(1,10), ones(1,5), zeros(1,10) ]';
nvox = length(signal); nsubj = 75; FWHM = 2;
store_thresh_rft = zeros(1, niters);
store_thresh_perm = zeros(1, niters);

rft_finds = 0;
perm_finds = 0;
for I = 1:niters
    modul(I,100);
    noisey_data = 10*randn([nvox,nsubj]) + signal;
    
    % RFT thresholding
    [im, store_thresh_rft(I)] = vRFT( noisey_data, FWHM );
    rft_finds = rft_finds + (sum(im(:)) > 0);
    
    % Permutation thresholidng
    store_thresh_perm(I) = perm_thresh(noisey_data, 'T', FWHM, NaN, 0);
    im_perm = mvtstat( noisey_data ) > store_thresh_perm(I);
    perm_finds = perm_finds + (sum(im_perm(:)) > 0);
end

rft_finds
perm_finds