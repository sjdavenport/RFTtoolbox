RSDmask = imgload('RSDmask_Beijing');
data = loadRSD(1:198, 'Cambridge', 4);

twoDdata = squeeze(data(:,:,50,:));
[bounds, bounded_mask] = mask_bounds(RSDmask);
twoDmask = squeeze(bounded_mask(:,:,50));

%%
tstat = mvtstat(data);
tstat_masked = tstat(bounded_mask);
histogram(tstat_masked);
max_tstat = max(tstat_masked(:))
tinv(max_tstat, 197)

%% 1D slice
D = 1; niters = 1000; data_slice = squeeze(twoDdata(:,30,:));
mask_slice = squeeze(twoDmask(:,30));
spfn = get_sample_fields( data_slice, mask_slice, D );
FWHM = 3; sample_size = 20; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
coverage1D = record_coverage( spfn, sample_size, params, niters)

%% 2D slice
D = 2; niters = 1000;
spfn = get_sample_fields( -twoDdata, twoDmask, D );
FWHM = 3; sample_size = 20; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
coverage2D = record_coverage( spfn, sample_size, params, niters)

%% Zeroed data
mean_zero_data = data - mean(data, 4);

twoDdata = squeeze(mean_zero_data(:,:,50,:));
[bounds, bounded_mask] = mask_bounds(RSDmask);
twoDmask = squeeze(bounded_mask(:,:,50));

D = 2; niters = 1000;
spfn = get_sample_fields( twoDdata, twoDmask, D );
FWHM = 3; sample_size = 20; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
coverage2D = record_coverage( spfn, sample_size, params, niters)

%% Adding scaling to ensure symmetry
twoDdata_sym = twoDdata;
twoDdata_sym(:,:,1:2:198) = -twoDdata(:,:,1:2:198);
tstat_2D = mvtstat(twoDdata_sym);
histogram(tstat_2D(logical(twoDmask)))

%% Scaled analysis
D = 2; niters = 1000;
spfn = get_sample_fields( twoDdata_sym, twoDmask, D );
FWHM = 3; sample_size = 40; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
record_coverage( spfn, sample_size, params, niters)

%% Paired analysis
paired_data = zeros([size(data, 1:3), 198/2]);
for I = 1:(198/2)
    paired_data(:,:,:,I) = data(:,:,:,2*I -1) - data(:,:,:,2*I);
end
twoDpaireddata = squeeze(paired_data(:,:,50,:));

D = 2; niters = 1000;
spfn = get_sample_fields( twoDpaireddata, twoDmask, D );
FWHM = 3; sample_size = 50; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
coverage2D = record_coverage( spfn, sample_size, params, niters)

%%
