%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the EEC functions
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D example
FWHM = 6; resadd = 1; nsubj = 50; nvox = 100;
lat_data = wnfield(nvox, nsubj);

[tcfield, cfields] = convfield_t_Field( lat_data, FWHM, resadd );
dcfields = convfield_Field( lat_data, FWHM, 1, resadd );
[L,L0] = LKC_voxmfd_est( cfields, dcfields );

thresholds = -5:0.1:5;
EEC_usingspm = EEC_calc( thresholds, L, L0, 'T', nsubj );
EEC_new = EEC( thresholds, L, L0, 'T', nsubj - 1 );

plot(thresholds, EEC_usingspm)
hold on
plot(thresholds, EEC_new, '--')

EEC_usingspm(10)
EEC_new(10)

%% %% 2D Examples
%% Simple 2D example
FWHM = 6; resadd = 1; nsubj = 50; Dim = [50,50];
lat_data = wnfield(Dim, nsubj);

[tcfield, cfields] = convfield_t_Field( lat_data, FWHM, resadd );
dcfields = convfield_Field( lat_data, FWHM, 1, resadd );
[L,L0] = LKC_voxmfd_est( cfields, dcfields );

thresholds = -5:0.1:5;

EEC_usingspm = EEC_calc( thresholds, L, L0, 'T', nsubj );
EEC_new = EEC( thresholds, L, L0, 'T', nsubj - 1 );

plot(thresholds, EEC_usingspm)
hold on
plot(thresholds, EEC_new, '--')

EEC_usingspm(10)
EEC_new(10)


%% %% 3D Examples
%% Simple 3D example
FWHM = 2; resadd = 1; nsubj = 50; Dim = [5,5,5];
lat_data = wnfield(Dim, nsubj);

[tcfield, cfields] = convfield_t_Field( lat_data, FWHM, resadd );
dcfields = convfield_Field( lat_data, FWHM, 1, resadd );
[L,L0] = LKC_voxmfd_est( cfields, dcfields );

thresholds = -5:0.1:5;

EEC_usingspm = EEC_calc( thresholds, L, L0, 'T', nsubj );
EEC_new = EEC( thresholds, L, L0, 'T', nsubj - 1 );

plot(thresholds, EEC_usingspm)
hold on
plot(thresholds, EEC_new, '--')

EEC_usingspm(10)
EEC_new(10)
