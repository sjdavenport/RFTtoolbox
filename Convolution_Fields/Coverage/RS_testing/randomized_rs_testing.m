data_path = 'C:\Users\12SDa\davenpor\Data\RestingStateData\EklundData\';
load([data_path 'UKB_2D_randomized'])

MNImask = imgload('MNImask');
MNImask_2D = MNImask(:,:,45);

%% Histogram of the data 
tstat_2D = mvtstat(im_store_2D(:,:,1:560));
histogram(tstat_2D(logical(MNImask_2D)))

%% Set Sample field function and parameters
D = 2; niters = 1000;
spfn = get_sample_fields( im_store_2D, MNImask_2D, D );
FWHM = 3; resadd = 1;

%% Record the coverage
sample_size = 200;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
coverage = record_coverage( spfn, sample_size, params, niters, 10)
ECcurveanal(coverage, MNImask_2D, sample_size, 0.1)

%% LKC calc
lat_data = Mask(spfn(100).lat_data); FWHM = 3; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );

[~,~,~,L] = vRFT(lat_data, params, 0);

cfield = convfield(lat_data,params);
dcfield = convfield(lat_data,params,1);

L_conv = LKC_voxmfd_est(cfield, dcfield)

HPE  = LKC_HP_est( cfield, 1000, 1);
HPEL = HPE.hatL'
            
%% Compare empirical EC curves to theory
sample_size = 20; lat_data = spfn(sample_size).lat_data;
tcfield = convfield_t(lat_data, params);
[empiricalECcurve, thresholds] = ECcurve( tcfield, [-6,6], 0.05);
L0 = EulerChar(MNImask_2D, 0.5, 2);
curve_conv = EEC( thresholds, L_conv, L0, 'T', sample_size-1 );
plot(thresholds, empiricalECcurve);
hold on
plot(thresholds, curve_conv)
legend('empirical', 'fit')

%% Average EC curves
sample_size = 50;
resadd = 5; FWHM = 3; params = ConvFieldParams( repmat(FWHM,1,D), resadd );
curve_store = computeECcurves( spfn, params, sample_size, niters );

%%
mean_EC_curve = mean(curve_store);
start_at = -6; initial_point = find(thresholds == start_at);
threshold_subset = thresholds(initial_point:end);
plot(threshold_subset, mean_EC_curve(initial_point:end) )
hold on
HPEL = HPE.hatL'
L0 = EulerChar(MNImask_2D, 0.5, 2)
curve_conv = EEC( threshold_subset, L_conv, L0, 'T', sample_size -1 );
plot(threshold_subset, curve_conv)
legend('average', 'Theory')

%% Compare to wnfield
spfn = @(nsubj) wnfield(logical(MNImask_2D), nsubj);
lat_data = spfn(100); FWHM = 3; resadd = 5;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
cfield = convfield(lat_data,params);
HPE  = LKC_HP_est( cfield, 1, 1);
HPEL = HPE.hatL'
% ECcurveanal(coverage, MNImask_2D, sample_size, 0.5, HPEL)

%% Compare empirical EC curve to theory
tcfield = convfield_t(lat_data, params);
[empiricalECcurve, thresholds] = ECcurve( tcfield, [-6,6], 0.05);
L0 = EulerChar(MNImask_2D, 0.5, 2);
curve_conv = EEC( thresholds, HPEL, L0, 'T', sample_size -1 );
plot(thresholds, empiricalECcurve);
hold on
plot(thresholds, curve_conv)
legend('empirical', 'fit')

%% View the data (Is the mask being applied correctly??)
lat_data = Mask(spfn(20).lat_data);
tstat = convfield_t(lat_data, params);
imagesc(tstat.field > -3)

%%
params = ConvFieldParams( repmat(FWHM,1,D), 0 );
lat_data = Mask(spfn(20).lat_data);
tstat = convfield_t(lat_data, params);
imagesc(tstat.*MNImask_2D)

%%
load([data_path 'UKB_3D_randomized'])
MNImask_3D = MNImask(35:45,40:50,35:45);

% Histogram of the data 
tstat_3D = mvtstat(im_store_3D);
histogram(tstat_3D(logical(MNImask_3D)))

%%
D = 3; niters = 1000;
spfn = get_sample_fields( im_store_3D, MNImask_3D, D );
FWHM = 3; sample_size = 20; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
coverage = record_coverage( spfn, sample_size, params, niters, 10)

ECcurveanal(coverage, MNImask_3D, sample_size, 0.1)
