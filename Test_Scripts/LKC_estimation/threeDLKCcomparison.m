%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script compares LKC estimation in 3D
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3D boundary sphere example - compare to theory
Dim = [5,5,5]; FWHM = 3; resadd = 1; D = length(Dim);
mask_sphere = true( Dim ); 
mask_sphere = bndry_voxels( mask_sphere, 'full');
pad = ceil(4*FWHM2sigma(FWHM));
mask_sphere_padded = logical( pad_vals( mask_sphere, pad) );
lat_data = wnfield(mask_sphere_padded, nsubj);

% Convolution L_1 (voxmndest)
cfield  = convfield_Field( lat_data, FWHM, 0, resadd );
dcfield = convfield_Field( lat_data, FWHM, 1, resadd );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );
L_conv

[~, detLambda] = Lambda_theory(FWHM, D);
sum(mask_sphere(:))*detLambda^(1/2)