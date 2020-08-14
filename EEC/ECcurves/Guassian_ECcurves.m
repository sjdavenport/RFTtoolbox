%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script compares the EC curves with the EEC curves for t-fields
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% 1D clusters
nvox = 100; niters = 1000; FWHM = 20;
D = 1; numberofclusters = 0; u = 1;
data = noisegen(nvox, niters, FWHM);
for I = 1:niters
    c_calc = clusterloc_1D(data(:,I)', u);
    numberofclusters = numberofclusters + c_calc.nC;
end
numberofclusters
average_number_of_clusters_above_u = numberofclusters/niters

% Compare to the theoretical prediction:
Gker_param = FWHM2sigma(FWHM);
Lambda_theory = Gker_param^(-2)/2;
L1 = nvox*Lambda_theory^(1/2); L0 = 1;
EEC_spm(u, L1, L0, 'Z')

%%
nvox = 100; niters = 1000; FWHM = 20;
D = 1; numberofclusters = 0; u = 1;
data = noisegen(nvox, 1, FWHM);
a = Field(size(data)); a.field = data;
[ curve, thresholds ] = ECcurve( a );
Gker_param = FWHM2sigma(FWHM);
Lambda_theory = Gker_param^(-2)/2;
L1 = nvox*Lambda_theory^(1/2); L0 = 1;
EEC = EEC_spm( thresholds, L1, L0, 'Z' );
plot(thresholds, curve)
hold on 
plot(thresholds, EEC)

%% Simple 1D example
nvox = 100; FWHM = 5; resadd = 1; nsubj = 200;
lat_data = wnfield(nvox, nsubj);

[tcfield, cfields] = convfield_t_Field( lat_data, FWHM, resadd );

dcfields = convfield_Field( lat_data, FWHM, 1, resadd );
[ curve, thresholds ] = ECcurve( tcfield );
[L,L0] = LKC_voxmfd_est( cfields, dcfields );
EEC_out = EEC( thresholds, L, L0, 'T', nsubj )
plot(thresholds, curve)
hold on 
plot(thresholds, EEC)

%% %% 2D Examples
%% Simple 2D example


%% %% 3D Examples
%% Simple 3D example