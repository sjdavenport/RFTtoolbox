nsubj_vec = 10:10:50;
tot_nsubj = 100;

nvox = 20; D = 1;

nsims = 100;

FWHM = 3; resadd = 1;  niters = 1000;
params = ConvFieldParams( FWHM, resadd );

coverage = struct();
for I = 1:nsims
    I
    for J = 1:length(nsubj_vec)
        J
        sample_size = nsubj_vec(J);
        spfn = @(nsubj) wtfield(nvox, nsubj, 3);
        coverage(I,J).coverage = record_coverage( spfn, sample_size, params, niters, 2);
    end
    save('./indep_coverage', 'coverage')
end
%%
nsubj_vec = 10:10:50;
tot_nsubj = 100;

nvox = 20; D = 1;

nsims = 100;

FWHM = 3; resadd = 1;  niters = 1000;
params = ConvFieldParams( FWHM, resadd );

coverage = struct();
for I = 1:nsims
    I
    for J = 1:length(nsubj_vec)
        J
        sample_size = nsubj_vec(J);
        spfn = @(nsubj) wnfield(nvox, nsubj);
        coverage(I,J).coverage = record_coverage( spfn, sample_size, params, niters, 2);
    end
    save('./indep_coverage_n', 'coverage')
end