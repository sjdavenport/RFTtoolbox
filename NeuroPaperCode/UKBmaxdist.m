function UKBmaxdist( FWHM, nsubj, RSfolder, niters  )
% UKBMAXDIST( FWHM, nsubj_vec, RSfolder, niters  ) obtain the distribution
% of the maximum of convolution fields
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  FWHM         the FWHM with which to smooth
%  nsubj        the number of subjects
%  RSfolder     a string giving the folder corresponding to the resting
%               state data
% Optional
%  niters   the number of iterations to do
%--------------------------------------------------------------------------
% OUTPUT
% Saves the subset, lattice maxima and convolution maxima in 
% [RFTboxloc,'NeuroPaperCode/Maxdist/md_FWHM_', num2str(FWHM), 'nsubj_', num2str(nsubj)]
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Make global the location of the RFTtoolbox
global RFTboxloc

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'niters', 'var' )
   % default option of opt1
   niters = 1000;
end

if ~exist('RSfolder', 'var')
    RSfolder = 'RS_2Block';
end

%%  Main Function Loop
%--------------------------------------------------------------------------
% Obtain the bounded mask
mask = imgload('MNImask');
[~,mask] = mask_bounds(mask);

% Use resadd = 0 for the initialization
resadd = 0;

% Initialize the storage vectors
store_sampleids = zeros(nsubj,niters);
npeaks = 8;
lat_maxima = zeros(npeaks, niters);
conv_maxima = zeros(npeaks, niters);

% Obtain the sample function that calls the data
spfn = get_sample_fields( RSfolder, mask );

% Obtain the params for the convolution
params = ConvFieldParams( [FWHM,FWHM,FWHM], resadd );

% Find the local maxima
for I = 1:niters
    sample_images = spfn(nsubj);
    store_sampleids(:,I) = sample_images.subset;
    
    tcfield = convfield_t( sample_images.lat_data, params );
    clear sample_images
    
    [ peak_est_locs, ~, peakvals ] = lmindices(tcfield.field, npeaks, tcfield.mask);
    lat_maxima(:,I) = peakvals';
    [~, max_tfield_at_lms] = findconvpeaks(sample_images.lat_data.field, FWHM, peak_est_locs, 'T', tcfield.mask);
    conv_maxima(:,I) = max_tfield_at_lms';
    save([RFTboxloc,'NeuroPaperCode/Maxdist/md_FWHM_', num2str(FWHM), 'nsubj_', num2str(nsubj)], 'lat_maxima', 'conv_maxima', 'store_sampleids')
end

end

