%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script compares the EC curves with the EEC curves for t-fields
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D example
nvox = 100; FWHM = 5; resadd = 1; nsubj = 20;
lat_data = wnfield(nvox, nsubj);

[tcfield, cfields] = convfield_t_Field( lat_data, FWHM, resadd );

dcfields = convfield_Field( lat_data, FWHM, 1, resadd );
[ curve, thresholds ] = ECcurve( tcfield, [-3,3]);
[L,L0] = LKC_voxmfd_est( cfields, dcfields );
EEC_val = EEC( thresholds, L, L0, 'T', nsubj )
plot(thresholds, curve)
hold on 
plot(thresholds, EEC_val)

%% %% 2D Examples
%% Simple 2D example
Dim = [100,100]; FWHM = 5; resadd = 1; nsubj = 20;
lat_data = wnfield(Dim, nsubj);
[tcfield, cfields] = convfield_t_Field( lat_data, FWHM, resadd );
dcfields = convfield_Field( lat_data, FWHM, 1, resadd );

[ curve, thresholds ] = ECcurve( tcfield, [-3,3]);
[L,L0] = LKC_voxmfd_est( cfields, dcfields );
EEC = EEC_spm( thresholds, L, L0, 'T', nsubj )
plot(thresholds, curve)
hold on 
plot(thresholds, EEC)

%% Note that SPM miscalcualtes the EC in 2D!! I.e:
[L,L0] = LKC_voxmfd_est( cfields, dcfields );
use_spm = 1;
EEC = EEC_spm( thresholds, L, L0, 'T', nsubj, use_spm )
plot(thresholds, curve)
hold on 
plot(thresholds, EEC)

%% %% 3D Examples
%% 3D example
Dim = [91,109,91]; FWHM = 3; resadd = 1; nsubj = 20;
lat_data = wnfield(Dim, nsubj); lat_data.mask = logical(imgload('MNImask'));
[tcfield, cfields] = convfield_t_Field( lat_data, FWHM, resadd );
dcfields = convfield_Field( lat_data, FWHM, 1, resadd );

[L,L0] = LKC_voxmfd_est( cfields, dcfields );

[ curve, thresholds ] = ECcurve( tcfield, [-5,5]);
EEC = EEC_calc( thresholds, L, L0, 'T', nsubj )
plot(thresholds, curve)
hold on 
plot(thresholds, EEC)

%% Masked Example
Dim = [91,109,91]; FWHM = 3; resadd = 1; nsubj = 20;
lat_data = wnfield(Dim, nsubj); lat_data.mask = logical(imgload('MNImask'));
[tcfield, cfields] = convfield_t_Field( lat_data, FWHM, resadd );
dcfields = convfield_Field( lat_data, FWHM, 1, resadd );

[L,L0] = LKC_voxmfd_est( cfields, dcfields );

[ curve, thresholds ] = ECcurve( tcfield, [-5,5]);
EEC = EEC_calc( thresholds, L, L0, 'T', nsubj )
plot(thresholds, curve)
hold on 
plot(thresholds, EEC)

%% Note that SPM also miscalculates the EC in 3D!! I.e:
[L,L0] = LKC_voxmfd_est( cfields, dcfields );
use_spm = 1;
EEC = EEC_spm( thresholds, L, L0, 'T', nsubj, use_spm )
plot(thresholds, curve)
hold on 
plot(thresholds, EEC)

%% SPM gives these dents!
EEC_using_spm = EEC_calc( thresholds, L, L0, 'T', nsubj, 1 )
EEC_nospm = EEC_calc( thresholds, L, L0, 'T', nsubj, 0 )
plot(thresholds, EEC_using_spm)
hold on 
plot(thresholds, EEC_nospm)
legend('SPM', 'Just multiply')