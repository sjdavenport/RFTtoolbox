%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script mastches L_D with theory
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1D
FWHM = 6; resadd = 1; nsubj = 50; nvox = 100;
lat_data = wnfield(nvox, nsubj);

% Convolution L_1 (voxmndest)
cfield  = convfield_Field( lat_data, FWHM, 0, resadd );
dcfield = convfield_Field( lat_data, FWHM, 1, resadd );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );
L_conv

% Theory
L_theory = LKC_isogauss_theory( FWHM, nvox )

% L1 via theory
[ ~, detLambda] = Lambda_theory( FWHM, 1);
top_L = nvox*detLambda^(1/2)

%% 2D sphere (with padding to ensure stationarity)
FWHM = 3; nsubj = 20; Dim = [5,5]; D = length(Dim);
resadd = 9; pad = ceil( 4 * FWHM2sigma( FWHM ) ); mask = true(Dim);
mask = bndry_voxels( mask, 'full');
padded_mask = logical( pad_vals( mask, pad) ); lat_data = wnfield(padded_mask, nsubj);

% Convolution L_1 (voxmndest) seems to work in this case
lat_masked = 0;
cfield  = convfield_Field( lat_data, FWHM, 0, resadd, lat_masked );
dcfield = convfield_Field( lat_data, FWHM, 1, resadd, lat_masked );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );
L_conv

% Theory
[ ~, detLambda] = Lambda_theory( FWHM, D);
LD_theory = sum(padded_mask(:))*detLambda^(1/2)

%% 3D boundary box example - compare to theory
Dim = [5,5,5]; FWHM = 3; resadd = 1; D = length(Dim);
mask = true( Dim ); pad = ceil(4*FWHM2sigma(FWHM));
mask_padded = logical( pad_vals( mask, pad) );
lat_data = wnfield(mask_padded, nsubj);

% Convolution L_1, note need to use lat_masked here to ensure that it
% smooths outside of the mask
lat_masked = 0;
cfield  = convfield_Field( lat_data, FWHM, 0, resadd, lat_masked );
dcfield = convfield_Field( lat_data, FWHM, 1, resadd, lat_masked );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );
L_conv

[ L_spm, L0 ] = LKC_SPM_est( FWHM, mask );
L_spm

% isogauss theory
LKC_theory = LKC_isogauss_theory( FWHM, Dim );
LKC_theory

[~, detLambda] = Lambda_theory(FWHM, D);
LD_theory = sum(mask_padded(:))*detLambda^(1/2)

%% 3D boundary sphere example - compare to theory
Dim = [5,5,5]; FWHM = 3; resadd = 1; D = length(Dim);
mask_sphere = true( Dim ); 
mask_sphere = bndry_voxels( mask_sphere, 'full');
pad = ceil(4*FWHM2sigma(FWHM));
mask_sphere_padded = logical( pad_vals( mask_sphere, pad) );
lat_data = wnfield(mask_sphere_padded, nsubj);

% Convolution L_1 (voxmndest)
lat_masked = 0;
cfield  = convfield_Field( lat_data, FWHM, 0, resadd, lat_masked );
dcfield = convfield_Field( lat_data, FWHM, 1, resadd, lat_masked );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );
L_conv

[~, detLambda] = Lambda_theory(FWHM, D);
sum(mask_sphere(:))*detLambda^(1/2)

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
