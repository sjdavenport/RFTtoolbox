data_path = 'C:\Users\12SDa\davenpor\Data\RestingStateData\EklundData\';
load([data_path 'UKB_2D'])
load([data_path, 'UKB_masks_2D'])
n_masks = size(mask_store, 3);

MNImask = imgload('MNImask');
MNImask_2D = MNImask(:,:,45);
new_mask = (1 - dilate_mask(mvtstat(im_store) < -3.6, 1)).*dilate_mask(MNImask_2D, -2);
imagesc(new_mask)

%% Adding scaling to ensure symmetry
scaled_im_store = im_store;
scaled_im_store(:,:,1:2:1331) = -im_store(:,:,1:2:1331);
tstat_2D = mvtstat(scaled_im_store);
histogram(tstat_2D(logical(new_mask)))

[bounds, bounded_mask] = mask_bounds(new_mask);
symdata = scaled_im_store(bounds{1}, bounds{2}, :);

%% Scaled analysis
D = 2; niters = 1000;
spfn = get_sample_fields( symdata, bounded_mask, D );
FWHM = 3; sample_size = 100; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
record_coverage( spfn, sample_size, params, niters)

%%
combined_mask = prod(mask_store(:,:,:) > 0.5, 3);
imagesc(combined_mask)

%%
tstat_2D = mvtstat(im_store);
histogram(tstat_2D(logical(MNImask_2D)))

%%
lat_data = Field(logical(new_mask));
lat_data.field = im_store;
params = ConvFieldParams( [3,3], 0);

tcfield = convfield_t( lat_data, params );
imagesc(Mask(tcfield))

nonzeros_entries = tcfield.field.*new_mask > 0 + (tcfield.field.*new_mask < 0)
histogram(tcfield.field(logical(new_mask)))

%%
imagesc(mvtstat(im_store) < -3.6)

%% Analysis with the new mask
D = 2; niters = 1000;
spfn = get_sample_fields( im_store, new_mask, D );
FWHM = 3; sample_size = 20; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
record_coverage( spfn, sample_size, params, niters)

%% Two-sample testing
tot_pairs = floor(1331/2);
new_data = zeros(91,109,tot_pairs-1);
for I = 1:tot_pairs-1
   new_data(:,:,I) =  im_store(:,:,2*I) - im_store(:,:,2*I+1);
end

[bounds, bounded_mask] = mask_bounds(new_mask);
new_data = new_data(bounds{1}, bounds{2}, :);

%% Analysing the difference data
im = new_data(:,:,13);
imagesc(mvtstat(smooth_data).*bounded_mask)

%% Histogram of data
lat_data = Field(bounded_mask);
lat_data.field = new_data;
FWHM = 3; resadd = 5;
params = ConvFieldParams( [FWHM, FWHM], resadd );
tcfield = convfield_t(lat_data, params);
% plot(squeeze(new_data(60,40,:)))
imagesc(tcfield.field.*tcfield.mask)
figure
histogram(tcfield.field(tcfield.mask))

%% Simulation histogram
lat_data = wnfield(bounded_mask, 600);
FWHM = 20; resadd = 0;
params = ConvFieldParams( [FWHM, FWHM], resadd );
tcfield = convfield_t(lat_data, params);
histogram(tcfield.field(bounded_mask))
% plot(squeeze(new_data(60,40,:)))

%%
tstat = mvtstat(new_data); tstat = tstat(logical(new_mask));
h = histogram(tstat);
plot(h.BinEdges(1:end-1) +h.BinWidth/2, h.Values/sum(new_mask(:))/h.BinWidth);
hold on
plot(-3:0.1:3, normpdf(-3:0.1:3, 0,1))

%% Analysis with the new mask
D = 2; niters = 1000;
spfn = get_sample_fields( -new_data, bounded_mask, D );
FWHM = 3; sample_size = 20; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
coverage = record_coverage( spfn, sample_size, params, niters)

%% Vsmall areas
D = 2; niters = 1000; data = new_data(40:50, 25:35, :); mask = bounded_mask(40:50,25:35);
spfn = get_sample_fields( data, mask, D );
FWHM = 3; sample_size = 100; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
record_coverage( spfn, sample_size, params, niters)

%% EC curve analysis
LKC_estimate = mean(coverage.storeLKCs,2)';
L0 = EulerChar(new_mask, 0.5,2);
[ curve, x ] = maxECcurve( coverage.convmaxima, 0.1 )
plot(x, curve);
curve_conv = EEC( x, LKC_estimate, L0, 'T', sample_size -1 );
hold on 
plot(x,curve_conv)
%% Simulation on the new mask
D = 2; niters = 1000;
data = normrnd(0,1,[size(bounded_mask), 600]);
spfn = get_sample_fields( data, bounded_mask, D );
FWHM = 3; sample_size = 200; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
record_coverage( spfn, sample_size, params, niters)

%% Analyzing some of the voxels
ts = squeeze(im_store(55,44,:));
plot(ts)
find(ts < -150)

mvtstat(ts')
mvtstat([ts(1:287); ts(289:end)]')

%%
imagesc(im_store(:,:,288))

%%
D = 2; niters = 1000;
spfn = get_sample_fields( im_store(:,:,1:n_masks), mask_store, D );
FWHM = 3; sample_size = 20; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
record_coverage( spfn, sample_size, params, niters)

%%
D = 1; niters = 1000;
spfn = get_sample_fields( -RS_1D_data, masks_1D, D );
FWHM = 3; sample_size = 20; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
coverage = record_coverage( spfn, sample_size, params, niters)

%% 1D LKC calculations
lat_data = spfn(20).lat_data;
cfield  = convfield_Field( lat_data, FWHM, 0, resadd );
LKC_HP_est( cfield, 1, 1 )

resadd = 11;
lat_data = spfn(20).lat_data;
cfield  = convfield_Field( lat_data, FWHM, 0, resadd );
LKC_HP_est( cfield, 1, 1 )

cfield  = convfield_Field( lat_data, FWHM, 0, resadd );
dcfield = convfield_Field( lat_data, FWHM, 1, resadd );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );

%% 1D analysis
%% Analysis with the new mask
D = 1; niters = 1000; data_slice = squeeze(new_data(:,30,:));
mask_slice = squeeze(bounded_mask(:,30));
spfn = get_sample_fields( data_slice, mask_slice, D );
FWHM = 3; sample_size = 20; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
coverage1D = record_coverage( spfn, sample_size, params, niters)

%% EC curve analysis
LKC_estimate = mean(coverage1D.storeLKCs,2)';
L0 = EulerChar(mask_slice, 0.5, D);
[ curve, x ] = maxECcurve( coverage1D.convmaxima, 0.1 )
plot(x, curve);
curve_conv = EEC( x, LKC_estimate, L0, 'T', sample_size -1 );
hold on 
plot(x,curve_conv)

%% 1D Gaussian analysis
D = 1; niters = 1000; data_slice = squeeze(new_data(:,30,:));
mask_slice = squeeze(bounded_mask(:,30));
spfn_orig = get_sample_fields( data_slice, mask_slice, D );
twenty_fields = spfn_orig(20).lat_data;
spfn = @(nsubj) Gmult(twenty_fields);
FWHM = 6; sample_size = 20; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
coverage1D = record_coverage( spfn, sample_size, params, niters)

%%
resadd = 11
FWHM = 6; sample_size = 20; 
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
plot(convfield_t(spfn(20), params))