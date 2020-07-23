%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests compares voxelwise permutation and RFT
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
fconv(lat_data.field, FWHM)
tcfield = convfield_t_Field( lat_data, FWHM, 0, 0 )

%% Simple 1D example (already with resadd = 1 RFT does better!)
signal = 2*[ zeros(1,10), ones(1,5), zeros(1,10) ]';
nvox = length(signal); nsubj = 75; FWHM = 1;
lat_data = 10*wnfield(nvox, nsubj) + signal;

% RFT thresholding
[im, threshold, ~, xvals] = vRFT( lat_data, FWHM );

% Permutation thresholidng
[im_perm, threshold_perm] = perm_thresh(lat_data.field, 'T', FWHM, NaN, 0);
% threshold_perm = perm_thresh(lat_data.field, 'T', FWHM, NaN, 1);
% im_perm = mvtstat( lat_data.field ) > threshold_perm;

% Plot results
subplot(3,1,1)
plot(xvals{1}, im); title('RFT Discoveries')
xlim(vec2lim(xvals{1}))
subplot(3,1,2)
plot(fconv(signal, FWHM)); title('Smoothed Signal');
xlim(vec2lim(xvals{1}))
subplot(3,1,3)
plot(im_perm); title('Permutation Discoveries')

threshold
threshold_perm
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
    lat_data = wnfield(nvox, nsubj);
    lat_data.field = 10*(lat_data.field) + signal;
    
    % RFT thresholding
    [im, store_thresh_rft(I)] = vRFT( lat_data, FWHM );
    rft_finds = rft_finds + (sum(im(:)) > 0);
    
    % Permutation thresholidng
    store_thresh_perm(I) = perm_thresh(lat_data.field, 'T', FWHM, NaN, 0);
    im_perm = mvtstat( lat_data.field ) > store_thresh_perm(I);
    perm_finds = perm_finds + (sum(im_perm(:)) > 0);
end

rft_finds
perm_finds

%% Simple 2D example
signal = zeros(11,11); signal(4:8,4:8) = 3;
Dim = size(signal); nsubj = 25; FWHM = 1;
lat_data = 10*wnfield(Dim, nsubj) + signal;

% RFT thresholding
[im, threshold, ~, xvals] = vRFT( lat_data, FWHM, 21 );

% Permutation thresholidng
[im_perm, threshold_perm] = perm_thresh(lat_data.field, 'T', FWHM, NaN, 0);

% Plot results
subplot(1,3,1)
imagesc(logical(im)); title('RFT Discoveries')
subplot(1,3,2)
imagesc(fconv(signal, FWHM));title('Smoothed Signal');
subplot(1,3,3)
imagesc(im_perm); title('Permutation Discoveries')


%% Higher noise example
signal = zeros(11,11); signal(4:8,4:8) = 3;
Dim = size(signal); nsubj = 25; FWHM = 3;
lat_data = 20*wnfield(Dim, nsubj) + signal;

% RFT thresholding
[im, threshold, ~, xvals] = vRFT( lat_data, FWHM, 21 );

% Permutation thresholidng
[im_perm, threshold_perm] = perm_thresh(lat_data.field, 'T', FWHM, NaN, 0);

% Plot results
subplot(1,3,1)
imagesc(logical(im)); title('RFT Discoveries')
subplot(1,3,2)
imagesc(fconv(signal, FWHM));title('Smoothed Signal');
subplot(1,3,3)
imagesc(im_perm); title('Permutation Discoveries')

%% multiple peaks
signal = zeros(11,11); signal(2:3,2:3) = 3;
signal(8:9,8:9) = 3; Dim = size(signal); nsubj = 25; FWHM = 3;
lat_data = 5*wnfield(Dim, nsubj) + signal;

% RFT thresholding
[im, threshold, ~, xvals] = vRFT( lat_data, FWHM, 21 );

% Permutation thresholidng
[im_perm, threshold_perm] = perm_thresh(lat_data.field, 'T', FWHM, NaN, 0);

% Plot results
subplot(1,3,1)
imagesc(logical(im)); title('RFT Discoveries')
subplot(1,3,2)
imagesc(fconv(signal, FWHM));title('Smoothed Signal');
subplot(1,3,3)
imagesc(im_perm); title('Permutation Discoveries')
