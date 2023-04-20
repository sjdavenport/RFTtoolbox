%% Gaussian random fields
nvox = 1000; niters = 1000; fwhm = 4; resadd = 0;
D = 1; numberofclusters = 0; u = 1.5; padding = round(4*FWHM2sigma(fwhm));
%Note that the padding is needed so that we have a stationary convolution
%field, since this is a stationary example. In general LKCs can be
%calculated under non-stationarity meaning that the edge effect doesn't
%matter.
lat_data = wfield(nvox + 2*padding, niters);
field_start_loc = spacep(padding+1, resadd)
field_end_loc = spacep(padding+nvox, resadd)

params = ConvFieldParams(fwhm, resadd);
smoothed_fields = convfield(lat_data, params);
vareverywhereest = var(smoothed_fields.field,0,2);
var_est = mean(vareverywhereest(field_start_loc:field_end_loc));
smoothed_fields = smoothed_fields*(1/sqrt(var_est)); % To ensure variance 1.

lm_count = 0;
for I = 1:niters
    modul(I,100)
    cfield = smoothed_fields.field(field_start_loc:field_end_loc,I)';
    c_calc = clusterloc_1D(cfield, u);
    numberofclusters = numberofclusters + c_calc.nC;
    
    local_maxima = lmindices(cfield, 'all');
    field_vals = cfield(local_maxima);
    lm_count = lm_count + sum(field_vals > u);
end
average_number_of_clusters_above_u = numberofclusters/niters
averagelms = lm_count/niters
Gker_param = FWHM2sigma(fwhm); Lambda_theory = Gker_param^(-2)/2;
L1 = nvox*Lambda_theory^(1/2); L0 = 1;
expectedEC = EEC( u, L1, L0 )

%%
c_0 = 1
c_1 = 1.5
EEC( c_1, L1, L0 )/EEC( c_0, L1, L0 )
exp(-(c_1^2-c_0^2)/2)

%% chi2 random fields
nvox = 1000; niters = 1000; fwhm = 4; resadd = 0;
D = 1; numberofclusters_0 = 0; numberofclusters_1 = 0;
u_1 = 1.5; u_0 = 1;
padding = round(4*FWHM2sigma(fwhm));
%Note that the padding is needed so that we have a stationary convolution
%field, since this is a stationary example. In general LKCs can be
%calculated under non-stationarity meaning that the edge effect doesn't
%matter.
lat_data = wfield(nvox + 2*padding, niters);
field_start_loc = spacep(padding+1, resadd)
field_end_loc = spacep(padding+nvox, resadd)

params = ConvFieldParams(fwhm, resadd);
smoothed_fields = convfield(lat_data, params);
vareverywhereest = var(smoothed_fields.field,0,2);
var_est = mean(vareverywhereest(field_start_loc:field_end_loc));
smoothed_fields = smoothed_fields*(1/sqrt(var_est)); % To ensure variance 1.

lm_count_0 = 0;
lm_count_1 = 0;
for I = 1:niters
    modul(I,100)
    cfield = (smoothed_fields.field(field_start_loc:field_end_loc,I)').^2;
    c_calc_1 = clusterloc_1D(cfield, u_1);
    c_calc_0 = clusterloc_1D(cfield, u_0);
    numberofclusters_0 = numberofclusters_0 + c_calc_0.nC;
    numberofclusters_1 = numberofclusters_1 + c_calc_1.nC;
    
    local_maxima = lmindices(cfield, 'all');
    field_vals = cfield(local_maxima);
    lm_count_1 = lm_count_1 + sum(field_vals > u_1);
    lm_count_0 = lm_count_0 + sum(field_vals > u_0);
end
averagelms_0 = lm_count_0/niters
average_number_of_clusters_above_u_0 = numberofclusters_0/niters

averagelms_1 = lm_count_1/niters
average_number_of_clusters_above_u_1 = numberofclusters_1/niters

%%
averagelms_1/averagelms_0
average_number_of_clusters_above_u_1/average_number_of_clusters_above_u_0
exp(-(u_1-u_0)/2)


%% chi2 random fields
nvox = 33; niters = 1000; fwhm = 2.76; fwhm = 6;
D = 1; numberofclusters_0 = 0; numberofclusters_1 = 0;
padding = round(4*FWHM2sigma(fwhm));
%Note that the padding is needed so that we have a stationary convolution
%field, since this is a stationary example. In general LKCs can be
%calculated under non-stationarity meaning that the edge effect doesn't
%matter.
field_start_loc_lat = spacep(padding+1, 0);
field_end_loc_lat = spacep(padding+nvox, 0);
resadd = 7;
field_start_loc_fine = spacep(padding+1, resadd);
field_end_loc_fine = spacep(padding+nvox, resadd);

lat_data = wfield(nvox + 2*padding, niters);

params = ConvFieldParams(fwhm, 0);
params_fine = ConvFieldParams(fwhm, resadd);

smoothed_fields = convfield(lat_data, params);
smoothed_fields_fine = convfield(lat_data, params_fine);

vareverywhereest = var(smoothed_fields.field,0,2);
var_est = mean(vareverywhereest(field_start_loc_lat:field_end_loc_lat));
smoothed_fields = smoothed_fields*(1/sqrt(var_est)); % To ensure variance 1.

vareverywhereest_fine = var(smoothed_fields_fine.field,0,2);
var_est_fine = mean(vareverywhereest_fine(field_start_loc_fine:field_end_loc_fine));
smoothed_fields_fine = smoothed_fields_fine*(1/sqrt(var_est_fine)); % To ensure variance 1.

vecofmaxima_lat = zeros(1,niters);
vecofmaxima_fine = zeros(1,niters);

for I = 1:niters
    loader( I, niters )
    lat_data = wfield(nvox + 2*padding, niters);

    cfield_lat = (smoothed_fields.field(field_start_loc_lat:field_end_loc_lat,I)').^2;
    cfield_fine = (smoothed_fields_fine.field(field_start_loc_fine:field_end_loc_fine,I)').^2;
    
    vecofmaxima_lat(I) = max(cfield_lat);
    vecofmaxima_fine(I) = max(cfield_fine);
end

mean(vecofmaxima_fine - vecofmaxima_lat)

%%
sigma_vals = 1:0.1:5;
lat_thresh = zeros(1, length(sigma_vals));
fine_thresh = zeros(1, length(sigma_vals));
for I = 1:length(sigma_vals)
    alpha_level = normcdf(-sigma_vals(I));
    lat_thresh(I) = prctile(vecofmaxima_lat, 100*(1-alpha_level));
    fine_thresh(I) = prctile(vecofmaxima_fine, 100*(1-alpha_level));
end
plot(sigma_vals, lat_thresh)
hold on
plot(sigma_vals, fine_thresh)
legend('lat thresh', 'fine thresh', 'Location','southeast')