data_path = 'C:\Users\12SDa\davenpor\Data\RestingStateData\EklundData\';
load([data_path 'UKB_2D'])

MNImask = imgload('MNImask');
MNImask_2D = MNImask(:,:,45);
new_mask = dilate_mask(MNImask_2D, -2);
imagesc(new_mask)

%% Histogram of the data 
tstat_2D = mvtstat(im_store);
histogram(tstat_2D(logical(new_mask)))

%% Two-sample testing
tot_pairs = floor(1331/2);
new_data = zeros(91,109,tot_pairs-1);
for I = 1:tot_pairs-1
   new_data(:,:,I) =  im_store(:,:,2*I) - im_store(:,:,2*I+1);
end

[bounds, bounded_mask] = mask_bounds(new_mask);
new_data = new_data(bounds{1}, bounds{2}, :);

%% Histogram of data
lat_data = Field(bounded_mask);
lat_data.field = new_data;
FWHM = 3; resadd = 1;
params = ConvFieldParams( [FWHM, FWHM], resadd );
tcfield = convfield_t(lat_data, params);
% plot(squeeze(new_data(60,40,:)))
imagesc(tcfield.field.*tcfield.mask)
figure
histogram(tcfield.field(tcfield.mask))

%% Record the coverage
D = 2; niters = 1000;
spfn = get_sample_fields( new_data, bounded_mask, D );
FWHM = 3; sample_size = 20; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
coverage = record_coverage( spfn, sample_size, params, niters)

%% EC curve analysis
LKC_estimate = mean(coverage.storeLKCs,2)';
L0 = EulerChar(new_mask, 0.5,2);
[ curve, x ] = maxECcurve( coverage.convmaxima, 0.1 )
plot(x, curve);
curve_conv = EEC( x, LKC_estimate, L0, 'T', sample_size -1 );
hold on 
plot(x,curve_conv)

%% Simulation on the new mask (just to make sure)
D = 2; niters = 1000;
data = normrnd(0,1,[size(bounded_mask), 600]);
spfn = get_sample_fields( data, bounded_mask, D );
FWHM = 3; sample_size = 200; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
record_coverage( spfn, sample_size, params, niters)

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

%% 1D Gaussian analysis (to see if non-Gaussianity is the issue, shocking have the same problem here too!)
D = 1; niters = 1000; data_slice = squeeze(new_data(:,30,:));
mask_slice = squeeze(bounded_mask(:,30));
spfn_orig = get_sample_fields( data_slice, mask_slice, D );
twenty_fields = spfn_orig(20).lat_data;
spfn = @(nsubj) Gmult(twenty_fields); %Gmult generates a Gaussian multiplier field
FWHM = 6; sample_size = 20; resadd = 1;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
coverage1D = record_coverage( spfn, sample_size, params, niters)
