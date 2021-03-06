%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script looks at the relationship between L_1 and the 0.05
%%%    threshold in 3D
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3D MNImask - compare to theory
FWHM = 3; resadd = 1; mask = logical(imgload('MNImask'));
pad = ceil(4*FWHM2sigma(FWHM)); nsubj = 20; D = 3;
mask_padded = logical( pad_vals( mask, pad) );
lat_data = wnfield(mask_padded, nsubj);

% Convolution L_1 (voxmndest)
lat_masked = 0;
cfield  = convfield_Field( lat_data, FWHM, 0, resadd, lat_masked );
dcfield = convfield_Field( lat_data, FWHM, 1, resadd, lat_masked );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );
L_conv

[~, detLambda] = Lambda_theory(FWHM, D);
sum(mask_padded(:))*detLambda^(1/2)

[ L_spm, L0 ] = LKC_SPM_est( FWHM, mask );

%%
L_1_vec = -500:25:0;
L = L_spm; sample_size = nsubj;

threshold = zeros(1, length(L_1_vec));
for I = 1:length(L_1_vec)
    L(1) = L_1_vec(I);
    resel_vec = LKC2resel(L, L0_conv);
    threshold(I) = spm_uc_RF_mod(0.05,[1,sample_size-1],'T',resel_vec,1);
end

plot(L_1_vec, threshold)

%% Threshold vs L0
mask = imgload('MNImask'); nsubj = 20; FWHM = 3;
[ L_spm, L0 ] = LKC_SPM_est( FWHM, mask );

L = L_spm; sample_size = nsubj;

resel_vec = LKC2resel(L, L0);
threshold = spm_uc_RF_mod(0.05,[1,sample_size-1],'T',resel_vec,1)

resel_vec = LKC2resel(L, L0+1);
threshold = spm_uc_RF_mod(0.05,[1,sample_size-1],'T',resel_vec,1)