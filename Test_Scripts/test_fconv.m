%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the fconv function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D example
lat_data = normrnd(0,1,1,100); FWHM = 3;
smoothed_fconv = fconv(lat_data, FWHM);
smoothed_spm = spm_conv(lat_data,FWHM);
plot(smoothed_spm); hold on; plot(smoothed_fconv)
legend('spm\_conv', 'fconv') 

%% 1D multiple subjects
nvox = 100; nsubj = 2; lat_data = normrnd(0,1,nvox,nsubj); FWHM = 3; D = 1;
smoothed_fconv = fconv(lat_data, FWHM, D)
smoothed_spm = zeros(nvox, nsubj);
for n = 1:nsubj
    smoothed_spm(:,n) = spm_conv(lat_data(:,n),FWHM);
end
plot(smoothed_spm, 'color',[0.85 0.325 0.0980]); hold on;
plot(smoothed_fconv, 'color',[0 0.447 0.7410]);

%% %% 2D Examples
%% Simple 2D example
lat_data = normrnd(0,1,25,25); FWHM = 3;
smoothed_fconv = fconv(lat_data, FWHM); 
smoothed_spm = spm_conv(lat_data,FWHM);
subplot(2,1,1)
surf(smoothed_fconv)
title('fconv')
subplot(2,1,2)
surf(smoothed_fconv)
title('SPM\_conv')

%% %% 3D Examples
Dim = [50,50,50]; lat_data = normrnd(0,1,Dim); halfDim = Dim(1)/2;
D = length(Dim); FWHM = 3; D = 3;
smoothed_spm = zeros(Dim);
spm_smooth(lat_data, smoothed_spm, FWHM);
smoothed_fconv = fconv(lat_data, FWHM);
sigma = FWHM2sigma(FWHM); truncation = ceil(6*sigma);
smoothed_fconv_spmkern = fconv(lat_data, @(x) spm_smoothkern(FWHM, x), D, truncation );
smoothed_cfield = convfield( lat_data, FWHM, 1, D, 0, 0);
plot(1:Dim(1),smoothed_fconv(:,halfDim,halfDim))
hold on 
plot(1:Dim(1),smoothed_spm(:,halfDim,halfDim))
plot(1:Dim(1),smoothed_cfield(:,halfDim,halfDim), '--')
plot(1:Dim(1),smoothed_fconv_spmkern(:,halfDim,halfDim), '--')
legend('fconv', 'SPM', 'convfield', 'fconv_smoothkern')

plot(-truncation:truncation, spm_smoothkern(FWHM, -truncation:truncation))
hold on
plot(-truncation:truncation, GkerMV(-truncation:truncation, FWHM))

% Compare speed to spm_smooth (much faster)
Dim = [50,50,50]; lat_data = normrnd(0,1,Dim);
tic; fconv(lat_data, FWHM); toc
tic; smoothed_spm = zeros(Dim);
tt = spm_smooth_mod(lat_data, smoothed_spm, FWHM); toc
