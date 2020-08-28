%% 1D
nvox = 30; Dim = [nvox, 1]; nsubj = 20; FWHM = 3;
mask = true(Dim);

[forman, keibel, conv] = compare_fwhm_ests( mask, nsubj, FWHM )

FWHM_vec = 1:0.5:6;
niters = 500;

forman.fwhm_ests = zeros(length(FWHM_vec), niters);
forman.Lambda_ests = zeros(length(FWHM_vec), niters);
kiebel.fwhm_ests = zeros(length(FWHM_vec), niters);
kiebel.Lambda_ests = zeros(length(FWHM_vec), niters);
conv.fwhm_ests = zeros(length(FWHM_vec), niters);
conv.Lambda_ests = zeros(length(FWHM_vec), niters);

for I = 1:length(FWHM_vec)
    FWHM = FWHM_vec(I);
    for J = 1:niters
        [forman_sim, keibel_sim, conv_sim] = compare_fwhm_ests( mask, nsubj, FWHM );
        forman.fwhm_ests(I, J) = forman_sim.fwhm;
        forman.Lambda_ests(I, J) = forman_sim.Lambda;
        forman.fwhm_ests(I, J) = forman_sim.fwhm;
        forman.Lambda_ests(I, J) = forman_sim.Lambda;
        forman.fwhm_ests(I, J) = forman_sim.fwhm;
        forman.Lambda_ests(I, J) = forman_sim.Lambda;
    end
end