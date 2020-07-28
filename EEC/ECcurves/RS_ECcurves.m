%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script calculates the EC curves of the Eklund resting state
%%%    data
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Beijing Data
RSDmask = imgload(['RSDmask_', 'Beijing']);
[ bounds, mask ] = mask_bounds( RSDmask ); nsubj = 50;
which_subs = randsample(198, nsubj, 0);
lat_data = Field(mask); lat_data.field = loadRSD(which_subs, 'Beijing');

%% Calcualte LKCs
FWHM = 2; resadd = 1; 

% Convolution L_1 (voxmndest)
[tcfield, cfields] = convfield_t_Field( lat_data, FWHM, resadd );
dcfield = convfield_Field( lat_data, FWHM, 1, resadd );
[L_conv,L0_conv] = LKC_voxmfd_est( cfields, dcfield )

% HPE (Worryingly L_3 can be negative here sometimes)
HPE  = LKC_HP_est( cfields, 1, 1);
newHPEestimate = HPE.hatL'
 
% SPM (Off of course due to the edge correction)
FWHM = est_smooth(cfields.field);
[ L_spm, L0 ] = LKC_SPM_est( FWHM, mask )

%% Compare EC curves (general i.e. for L1 examination)
[ curve, thresholds ] = ECcurve( tcfield, [-5,5], 0.001);

% EEC_spm = EEC_calc( thresholds, L_spm, L0, 'T', nsubj ) SPM is very wrong!
EEC_hpe = EEC_calc( thresholds, newHPEestimate, L0, 'T', nsubj )
EEC_conv = EEC_calc( thresholds, L_conv, L0, 'T', nsubj )
EEC_mix = EEC_calc( thresholds, [newHPEestimate(1), L_conv(2), L_conv(3)], L0, 'T', nsubj )

plot(thresholds, curve)
hold on 
plot(thresholds, EEC_hpe)
plot(thresholds, EEC_conv)
plot(thresholds, EEC_mix)
legend('Observed', 'HPE', 'Conv', 'mix')

%% Compare EC curves (general i.e. for L1 examination)
[ curve, thresholds ] = ECcurve( tcfield, [3,6], 0.01);

EEC_spm = EEC_calc( thresholds, L_spm, L0, 'T', nsubj )
EEC_hpe = EEC_calc( thresholds, newHPEestimate, L0, 'T', nsubj )
EEC_conv = EEC_calc( thresholds, L_conv, L0, 'T', nsubj )
EEC_mix = EEC_calc( thresholds, [newHPEestimate(1), L_conv(2), L_conv(3)], L0, 'T', nsubj )

plot(thresholds, curve)
hold on
plot(thresholds, EEC_spm)
hold on 
plot(thresholds, EEC_hpe)
plot(thresholds, EEC_conv)
plot(thresholds, EEC_mix)
legend('Observed', 'SPM', 'HPE', 'Conv')

%% Average EC curve
niters = 500;
RSDmask = imgload(['RSDmask_', 'Beijing']);
[ bounds, mask ] = mask_bounds( RSDmask ); nsubj = 50;
which_subs = randsample(198, nsubj, 0);
lat_data = Field(mask); lat_data.field = loadRSD(which_subs, 'Beijing');
tic
[tcfield, cfields] = convfield_t_Field( lat_data, FWHM, resadd );
toc


%%


