data_path = 'C:\Users\12SDa\davenpor\Data\RestingStateData\EklundData\';
load([data_path 'UKB_2D_randomized'])

MNImask = imgload('MNImask');
MNImask_2D = MNImask(:,:,45);

%% Histogram of the data 
tstat_2D = mvtstat(im_store_2D);
histogram(tstat_2D(logical(MNImask_2D)))

%% Record the coverage
D = 2; niters = 1000;
spfn = get_sample_fields( im_store_2D, MNImask_2D, D );

nsubj_vec = 10:20:110;
store_coverage = struct();
FWHM = 6; resadd = 1;
params = ConvFieldParams([FWHM, FWHM], resadd);

for I = 1:length(nsubj_vec)
    nsubj = nsubj_vec(I);
    store_coverage.(['nsubj_',num2str(nsubj)]) = record_coverage(spfn, nsubj, params, niters);
    save(['./store_coverage_odd_FWHM', num2str(FWHM)], 'store_coverage')
end

%% Record the coverage
D = 2; niters = 1000;
spfn = get_sample_fields( im_store_2D, MNImask_2D, D );

nsubj_vec = 20:20:110;
store_coverage = struct();
FWHM = 3; resadd = 1;
params = ConvFieldParams([FWHM, FWHM], resadd);

for I = 1:length(nsubj_vec)
    nsubj = nsubj_vec(I);
    store_coverage.(['nsubj_',num2str(nsubj)]) = record_coverage(spfn, nsubj, params, niters);
    save('./store_coverage_even', 'store_coverage')
end

%% EC curve analysis
lat_data = Field(logical(MNImask_2D));
lat_data.field = im_store_2D;
FWHM = 3; resadd = 1;
params = ConvFieldParams([FWHM,FWHM], resadd);

cfield = convfield(lat_data, params);
dcfield = convfield(lat_data, params, 1);
L_conv = LKC_voxmfd_est(cfield, dcfield)

L_HPE = LKC_HP_est( cfield, 1000, 1)

%%
sample_size = 110;
coverage = store_coverage.(['nsubj_',num2str(sample_size)]);
ECcurveanal(coverage, MNImask_2D, sample_size, 0.1, L_conv)
ECcurveanal(coverage, MNImask_2D, sample_size, 0.1, L_HPE.hatL')

%%
global RFTboxloc
load([RFTboxloc, 'Convolution_Fields/Coverage/RS_testing/store_coverage_odd.mat'])
sample_size = 110;
coverage = store_coverage.(['nsubj_',num2str(sample_size)]);
ECcurveanal( coverage, MNImask_2D, sample_size, 0.1)