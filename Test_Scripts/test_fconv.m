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
figure, clf,
plot(smoothed_spm); hold on; plot(smoothed_fconv, '--')
legend('spm\_conv', 'fconv') 

%% 1D multiple subjects
nvox = 100; nsubj = 2; FWHM = 3; D = 1;
lat_data = normrnd( 0, 1, nvox, nsubj );
smoothed_fconv = fconv( lat_data, FWHM, D );
smoothed_spm = zeros( nvox, nsubj );
for n = 1:nsubj
    smoothed_spm(:,n) = spm_conv( lat_data(:,n), FWHM );
end
plot( smoothed_spm, 'color',  [ 0 0.447 0.7410 ] ); hold on;
plot( smoothed_fconv, '--', 'color', [ 0.85 0.325 0.0980 ] );

%% %% 2D Examples
%% % Simple 2D example using numeric FWHM input
lat_data = normrnd( 0, 1, 25, 25 ); FWHM = 5;
smoothed_fconv = fconv( lat_data, FWHM ); 
smoothed_spm = spm_conv( lat_data, FWHM );
subplot(2,1,1)
surf(smoothed_fconv)
title('fconv')
subplot(2,1,2)
surf(smoothed_fconv)
title('SPM\_conv')

%% % Simple 2D example using numeric vector FWHM input
% Generate lattice data
lat_data = normrnd( 0, 1, 25, 25 );
% FWHM  for smoothing in each direction
FWHM = [3, 7];
% change the frequency where the kernel is evaluated. Voxels are considered
% dx units away
dx = 0.5;
% Smooth the data
smoothed_fconv = fconv( lat_data, FWHM ); 
figure, clf,
surf(smoothed_fconv)
title('fconv defaults')

%% % Same result but filling the arguments manually
% Create a Sep_Kernel objcet representing the Gaussian kernel
K = SepKernel( 2, FWHM )
% Smooth the data
smoothed_fconv = fconv( lat_data, K.kernel, 2, K.truncation, 1, K.adjust );
figure, clf,
surf(smoothed_fconv)
title('fconv manually filled')

%% % Change of size of voxels
% change the frequency where the kernel is evaluated. Voxels are considered
% dx units away
dx = 0.5;
% Smooth the data
smoothed_fconv = fconv( lat_data, K.kernel, 2, K.truncation, dx, K.adjust );
figure, clf,
surf(smoothed_fconv)
title('fconv voxels changed to one half')

%% %% 3D Examples
Dim = [50,50,50]; lat_data = normrnd(0,1,Dim); halfDim = Dim(1)/2;
FWHM = 3; D = 3;
smoothed_spm = zeros(Dim);
spm_smooth(lat_data, smoothed_spm, FWHM);
smoothed_fconv = fconv(lat_data, FWHM);
sigma = FWHM2sigma(FWHM); truncation = ceil(6*sigma);
smoothed_fconv_spmkern = fconv(lat_data, @(x) spm_smoothkern(FWHM, x), D, truncation );
smoothed_cfield = convfield( lat_data, FWHM, 0, D );
figure, clf,
plot(1:Dim(1),smoothed_fconv(:,halfDim,halfDim))
hold on 
plot(1:Dim(1),smoothed_spm(:,halfDim,halfDim))
plot(1:Dim(1),smoothed_cfield(:,halfDim,halfDim), '--')
plot(1:Dim(1),smoothed_fconv_spmkern(:,halfDim,halfDim), '--')
legend('fconv', 'SPM', 'convfield', 'fconv_smoothkern')

figure, clf,
plot(-truncation:truncation, spm_smoothkern(FWHM, -truncation:truncation))
hold on
plot(-truncation:truncation, GkerMV(-truncation:truncation, FWHM))
title("Difference in Gaussian kernel and spm12 kernel")

% Compare speed to spm_smooth (much faster)
Dim = [50,50,50]; lat_data = normrnd(0,1,Dim);
tic; fconv(lat_data, FWHM); toc
tic; smoothed_spm = zeros(Dim);
tt = spm_smooth_mod(lat_data, smoothed_spm, FWHM); toc
