%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the reimannian metric estimation to see why 
%%%    there's not much difference to SPM in the simple example considered
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot estimates of the reimannian metric
FWHM = 3; sample_size = 50; Dim = [5,5]; mask = true(Dim); resadd = 11;
nsubj = 50; lat_data = normrnd(0,1,[Dim, nsubj]);

ghat = Lambda_numeric_est( lat_data, FWHM, resadd );
first_entry = ghat(:,:,1,1);
mean(ghat(:))
subplot(2,1,1);imagesc(first_entry); title('Lambda num est')
ghat = Lambda_conv_est( lat_data, FWHM, resadd );
first_entry = ghat(:,:,1,1);
mean(ghat(:))
subplot(2,1,2); imagesc(first_entry); title('Lambda conv est')

%% Compare stationary estimate to stationary theory (off course there is an edge effect so not suprising that there is a difference)
sigma = FWHM2sigma(FWHM); D = length(Dim);
Lambda_theory = diag(repmat(sigma^(-2),1,D))/2

smooth_data = fconv(lat_data, FWHM);
fwhm_spm_est = est_smooth(smooth_data, mask);
sigma_spm_est = FWHM2sigma(fwhm_spm_est);
Lambda_spm_est = diag(sigma_spm_est.^(-2))/2

%% Comparing different LKC calculations
FWHM = 3; sample_size = 50; Dim = [5,5]; mask = true(Dim); resadd = 11;
nsubj = 50; lat_data = normrnd(0,1,[Dim, nsubj]);

lkcs_conv = LKC_conv_est( lat_data, mask, FWHM, resadd );
lkcs_conv.hatL
lkcs_conv.L0

smooth_data = fconv(lat_data, FWHM);
lkcs_hpe = LKC_HP_est( smooth_data, mask, 1 );
lkcs_hpe.hatL
lkcs_hpe.L0

fwhm_spm_est = est_smooth(smooth_data, mask);
resels_est = spm_resels(fwhm_spm_est,Dim, 'B');
lkcs_spm_est = resel2LKC(resels_est)

spm_resels(FWHM,Dim, 'B');
resels = spm_resels(FWHM,Dim, 'B');
lkcs_spm_truefwhm = resel2LKC(resels)

%% Comparing thresholds (SPM gives a slightly higher threshold)
FWHM = 3; sample_size = 50; Dim = [5,5]; mask = true(Dim); resadd = 11;
nsubj = 50; lat_data = normrnd(0,1,[Dim, nsubj]);

lkcs_conv = LKC_conv_est( lat_data, mask, FWHM, resadd );
resels_conv = LKC2resel(lkcs_conv);
threshold = spm_uc_RF_mod(0.05,[1,sample_size-1],'T',resels_conv,1)

smooth_data = fconv(lat_data, FWHM);
fwhm_spm_est = est_smooth(smooth_data, mask);
resels_est = spm_resels(fwhm_spm_est,Dim, 'B');
threshold_spm = spm_uc_RF_mod(0.05,[1,sample_size-1],'T',resels_est,1)
