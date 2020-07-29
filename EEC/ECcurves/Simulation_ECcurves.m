%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script generates 3D simulation EC curves on brain masks
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% EC curve (Gaussian sim)
FWHM = 10; nsubj = 20; resadd = 1; lat_data = wnfield(RSDmask, nsubj);

% Convolution L_1 (voxmndest)
[tcfield, cfields] = convfield_t_Field( lat_data, FWHM, resadd );
dcfield = convfield_Field( lat_data, FWHM, 1, resadd );
[L_conv,L0_conv] = LKC_voxmfd_est( cfields, dcfield )

% HPE (Worryingly L_3 can be negative here sometimes)
HPE  = LKC_HP_est( cfields, 1, 1);
L_hpe = HPE.hatL'
 
% SPM (Off of course due to the edge correction)
[ L_spm, L0 ] = LKC_SPM_est( FWHM, mask )

[ curve, thresholds ] = ECcurve( tcfield, [3,6], 0.01);

EEC_spm = EEC_calc( thresholds, L_spm, L0, 'T', nsubj )
EEC_hpe = EEC_calc( thresholds, L_hpe, L0, 'T', nsubj )
EEC_conv = EEC_calc( thresholds, L_conv, L0, 'T', nsubj )
EEC_mix = EEC_calc( thresholds, [L_hpe(1), L_conv(2), L_conv(3)], L0, 'T', nsubj )

plot(thresholds, curve)
hold on
plot(thresholds, EEC_spm)
hold on 
plot(thresholds, EEC_hpe)
plot(thresholds, EEC_conv)
plot(thresholds, EEC_mix)
legend('Observed', 'SPM', 'HPE', 'Conv', 'Mix')

%% EC curve (Gaussian sim, just SPM)
FWHM = 10; nsubj = 20; resadd = 1; lat_data = wnfield(RSDmask, nsubj);

% Convolution L_1 (voxmndest)
tcfield = convfield_t_Field( lat_data, FWHM, resadd );

% SPM (Off of course due to the edge correction)
[ L_spm, L0 ] = LKC_SPM_est( FWHM, mask );

[ curve, thresholds ] = ECcurve( tcfield, [-5,5], 0.01);

EEC_spm = EEC_calc( thresholds, L_spm, L0, 'T', nsubj );

plot(thresholds, curve)
hold on
plot(thresholds, EEC_spm)
legend('Observed', 'SPM')

