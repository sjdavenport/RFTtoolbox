%% 1D coverage
FWHM_vec = 1:0.5:6;
nsubj_vec = 10:10:100;

Dim = 100;
mask = true([100,1]);

store_coverage( Dim, mask, FWHM_vec, nsubj_vec, 1)

%% 1D with mask

%% 3D coverage (small scale)
FWHM_vec = 3:6;
nsubj_vec = 25;

Dim = [5,5,5];
mask = true(Dim);

resadd = 3;

store_coverage( Dim, mask, FWHM_vec, nsubj_vec, 0, resadd)

%% 3D - test
FWHM_vec = 1:6; nsubj_vec = 25; Dim = [5,5,5];
mask = true(Dim); mask(2:4,2:4,2:4) = 0;
resadd = 3; niters = 1000;

store_coverage( Dim, mask, FWHM_vec, nsubj_vec, 0, resadd, niters, '3D_small_mask')


%% 3D - small - with mask
FWHM_vec = 1:6; nsubj_vec = 25; Dim = [5,5,5];
mask = true(Dim); mask(2:4,2:4,2:4) = 0;
resadd = 3; niters = 1000;

store_coverage( Dim, mask, FWHM_vec, nsubj_vec, 0, resadd, niters, '3D_small_mask')

%% 3D - small - with mask - highres
FWHM_vec = 1:6; nsubj_vec = 25; Dim = [5,5,5];
mask = true(Dim); mask(2:4,2:4,2:4) = 0;
resadd = 9; niters = 1000;

store_coverage( Dim, mask, FWHM_vec, nsubj_vec, 0, resadd, niters, '3D_small_mask_highres')

%% 3D - small - with mask - highres
FWHM_vec = 3; nsubj_vec = 25; Dim = [91,109,91];
mask = logical(imgload('MNImask'));
resadd = 1; niters = 1;

tic
store_coverage( Dim, mask, FWHM_vec, nsubj_vec, 0, resadd, niters, '3D_MNImask_example')
toc