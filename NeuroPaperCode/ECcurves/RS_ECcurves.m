%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script calculates the EC curves of the resting state data
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dim = [91,109,91]; FWHM = 3; resadd = 1; nsubj = 20;
lat_data = wnfield(Dim, nsubj);
[tcfield, cfields] = convfield_t_Field( lat_data, FWHM, resadd );
dcfields = convfield_Field( lat_data, FWHM, 1, resadd );

[ curve, thresholds ] = ECcurve( tcfield, [-5,5]);
[L,L0] = LKC_voxmfd_est( cfields, dcfields );
EEC = EEC_calc( thresholds, L, L0, 'T', nsubj )
plot(thresholds, curve)
hold on 
plot(thresholds, EEC)

