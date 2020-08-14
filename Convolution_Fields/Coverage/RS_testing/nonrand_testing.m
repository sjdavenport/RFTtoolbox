data_path = 'C:\Users\12SDa\davenpor\Data\RestingStateData\EklundData\';
load([data_path 'UKB_2D'])

MNImask = imgload('MNImask');
MNImask_2D = MNImask(:,:,45);

data_mean = mean(im_store, 3);
im_store = im_store - data_mean;

%% Histogram of the demeaned data 
tstat_2D = mvtstat(im_store(:,:,1:200));
histogram(tstat_2D(logical(MNImask_2D)))

%% Record the coverage
D = 2; niters = 1000;
spfn = get_sample_fields( im_store, MNImask_2D, D );

nsubj_vec = 50:50:500;
store_coverage = struct();
FWHM = 3; resadd = 1;
params = ConvFieldParams([FWHM, FWHM], resadd);

for I = 1:length(nsubj_vec)
    nsubj = nsubj_vec(I)
    store_coverage.(['nsubj_',num2str(nsubj)]) = record_coverage(spfn, nsubj, params, niters);
    save(['./store_coverage_nr_FWHM', num2str(FWHM)], 'store_coverage')
end

%%
global RFTboxloc
load([RFTboxloc, 'Convolution_Fields/Coverage/RS_testing/store_coverage_nr_FWHM3.mat'])
sample_size = 250;
coverage = store_coverage.(['nsubj_',num2str(sample_size)]);
ECcurveanal( coverage, MNImask_2D, sample_size, 0.1)