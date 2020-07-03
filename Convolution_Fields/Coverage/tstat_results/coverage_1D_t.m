%% 1D coverage
FWHM_vec = 1:0.5:6;
nsubj_vec = 10:10:100;

Dim = 100;
mask = true([100,1]);

store_coverage( Dim, mask, FWHM_vec, nsubj_vec, 1)

%% 3D coverage (small scale)
FWHM_vec = 3:6;
nsubj_vec = 25;

Dim = [5,5,5];
mask = true(Dim);

resadd = 3;

store_coverage( Dim, mask, FWHM_vec, nsubj_vec, 0, resadd)