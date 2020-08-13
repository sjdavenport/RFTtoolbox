data_path = 'C:\Users\12SDa\davenpor\Data\RestingStateData\EklundData\';
load([data_path 'UKB_2D_randomized'])

MNImask = imgload('MNImask');
MNImask_2D = MNImask(:,:,45);
D = 2; niters = 500;

%% Multiplier fields
base_fields = Field(logical(MNImask_2D));
base_fields.field = im_store_2D(:,:,21:40);
spfn = @(nsubj)  Mask( multiplier_field( base_fields, nsubj ) );

%% Gaussian analysis (Multiplier field)
coverage2D = record_coverage( spfn, sample_size, params, niters)
ECcurveanal( coverage2D, MNImask_2D, sample_size )

%% Generate multiplier sample for bootstrapping
rng(300)
base_fields = Field();
base_fields.field = im_store_2D(:,:,21:40);
spfn = @(nsubj) multiplier_field( base_fields, nsubj );
sample_data = spfn(500).field;

%% Multiplier fields
spfn = get_sample_fields( sample_data, MNImask_2D, D );
sample_size = 20;
FWHM = 3; resadd = 1;
params = ConvFieldParams([FWHM,FWHM], resadd); niters = 1000;
coverage2D = record_coverage( spfn, sample_size, params, niters, 10)
ECcurveanal( coverage2D, MNImask_2D, sample_size )

%% Multiplier testing with extra mask
sample_data_zeroed = sample_data.*MNImask_2D;
spfn = get_sample_fields( sample_data_zeroed, MNImask_2D, D );
sample_size = 20;
FWHM = 3; resadd = 1;
params = ConvFieldParams([FWHM,FWHM], resadd); niters  = 1000;
coverage2D = record_coverage( spfn, sample_size, params, niters)
ECcurveanal( coverage2D, MNImask_2D, sample_size )

%% LKC estimates
lat_data = Mask(spfn(100)); FWHM = 3; resadd = 5;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
cfield = convfield(lat_data,params);
HPE  = LKC_HP_est( cfield, 1, 1);
HPEL = HPE.hatL'

%% Compare empirical EC curve to theory
sample_size = 20; lat_data = spfn(sample_size);
tcfield = convfield_t(lat_data, params);
[empiricalECcurve, thresholds] = ECcurve( tcfield, [-6,6], 0.05);
L0 = EulerChar(MNImask_2D, 0.5, 2);
curve_conv = EEC( thresholds, HPEL, L0, 'T', sample_size -1 );
plot(thresholds, empiricalECcurve);
hold on
plot(thresholds, curve_conv)
legend('empirical', 'fit')

%% Average EC curves
resadd = 5; FWHM = 3; params = ConvFieldParams( repmat(FWHM,1,D), resadd );
[curve_store, thresholds] = computeECcurves( spfn, params, sample_size, niters );

%% Plot average EC curve vs EEC estimate
plot(thresholds, mean(curve_store))
hold on
HPEL = HPE.hatL'
L0 = EulerChar(MNImask_2D, 0.5, 2)
curve_conv = EEC( thresholds, HPEL, L0, 'T', sample_size -1 );
plot(thresholds, curve_conv)
legend('average', 'HPE')
 