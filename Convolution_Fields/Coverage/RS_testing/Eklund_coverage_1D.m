RSDmask = imgload('RSDmask_Beijing');
data = loadCF(1:198, 'E4', 6);

twoDdata = squeeze(data(:,:,50,:));
MNImask = imgload('MNImask');
[bounds, bounded_mask] = mask_bounds(RSDmask);
twoDmask = squeeze(bounded_mask(:,:,50));

%% Set Sample field function and parameters
D = 2; niters = 1000;
spfn = get_sample_fields( -twoDdata, twoDmask, D );
FWHM = 3; resadd = 1;

%% Record the coverage
sample_size = 40;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
coverage = record_coverage( spfn, sample_size, params, niters )
ECcurveanal(coverage, twoDmask, sample_size, 0.1)