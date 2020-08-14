%%
load([data_path 'UKB_3D_randomized'])
MNImask_3D = MNImask(35:45,40:50,35:45);

% Histogram of the data 
tstat_3D = mvtstat(im_store_3D);
histogram(tstat_3D(logical(MNImask_3D)))

%%
D = 3; niters = 1000;
spfn = get_sample_fields( im_store_3D, MNImask_3D, D );
FWHM = 3; sample_size = 50; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
coverage = record_coverage( spfn, sample_size, params, niters, 10)
ECcurveanal(coverage, MNImask_3D, sample_size, 0.1)
