%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the EEC functions
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D example
FWHM = 6; resadd = 1; nsubj = 10; nvox = 100;
lat_data = wfield(nvox, nsubj);
params = ConvFieldParams(FWHM, resadd);

[tcfield, cfields] = convfield_t(lat_data, params);
dcfields = convfield(lat_data, params, 1);
[L,L0] = LKC_voxmfd_est( cfields, dcfields );

thresholds = -5:0.1:5;
EEC_usingspm = EEC_calc( thresholds, L, L0, 'T', nsubj );
EEC_new = EEC( thresholds, L, L0, 'T', nsubj - 1 );

plot(thresholds, EEC_usingspm)
hold on
plot(thresholds, EEC_new, '--')

EEC_usingspm(10)
EEC_new(10)

%% Expected number of maxima
niters = 1000;
FWHM = 7;
resadd = 5;
nvox = 100;
params = ConvFieldParams(FWHM, resadd);
n_maxima = 0;
threshold = 4;
for I = 1:niters
    modul(I, 100)
    lat_data = wfield(nvox, nsubj);
    f = convfield_t(lat_data, params);
    n_maxima = n_maxima + length(lmindices(f.field, 'all', f.field > threshold));
end
n_maxima/niters
bernstd(n_maxima/niters, niters)
EEC_calc( threshold, L, L0, 'T', nsubj )

%% 1D stationary example
FWHM = 6; resadd = 1; nsubj = 100; nvox = 100;
params = ConvFieldParams(FWHM, resadd);
trunc = ceil(4*FWHM2sigma(FWHM));
noise = wfield(nvox + 2*trunc, nsubj);
smooth_f = convfield(noise, params);
start = floor((trunc/(nvox+2*trunc))*smooth_f.masksize(1));
stat_f = smooth_f(start:smooth_f.masksize(1)-start);
smooth_f_deriv = convfield(noise, params, 1);
stat_f_derivs = smooth_f_deriv(start:smooth_f.masksize(1) -start);

% compute the tstat but no additional smoothing
[ L, L0, Lambda ] = LKC_stationary_est( stat_f, stat_f_derivs )
FWHM2Lambda(FWHM,1)

thresholds = -5:0.1:5;
EEC_usingspm = EEC_calc( thresholds, L, L0, 'T', nsubj );
EEC_new = EEC( thresholds, L, L0, 'T', nsubj - 1 );

plot(thresholds, EEC_usingspm)
hold on
plot(thresholds, EEC_new, '--')

EEC_usingspm(10)
EEC_new(10)

%% Expected number of maxima


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
