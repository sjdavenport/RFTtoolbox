function coverage = record_coverage( spfn, sample_size, Kernel, resadd, mask, niters, lkc_est_version, do_spm )
% RECORD_COVERAGE( data, FWHM, mask, B, sample_size ) estimates the coverage
% provided by a variety of RFT implementations including non-stationary and
% stationary convolution and lattice versions.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  spnfn         a function handle
%  sample_size   the size of each sample to be sampled from the data
%  FWHM          the applied FWHM of the Gaussian Kernel in each direction
%          (we smooth with an istropic Kernel as is commonly done in practice)
% Optional
%  resadd       a non-negative integer giving the resolution increase.
%               Default is 1.
%  mask          a 0/1 array of size Dim which provides a mask of the data
%               the default to use is no mask i.e. 1
%  niters        the number of resamples of the data to do
%  lkc_est_version      either 'conv' or 'hpe'. Default is 'conv'
%  do_spm       additionally calculate the lkcs using SPM (i.e. under
%               stationarity
%--------------------------------------------------------------------------
% OUTPUT
%
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------


%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
sample_image_size = size(spfn(1));
% Dim = sample_image_size(1:end-1); %This returns the image dimension by construction.
if sample_image_size(1) == 1
    D = 1;
    Dim = sample_image_size(2);
elseif sample_image_size(2) == 1
    D = 1;
    Dim = sample_image_size(1);
else
    Dim = sample_image_size;
    D = length(Dim);
end

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'niters', 'var' )
    % default option of niters
    niters = 1000;
end

% Set the default mask value
if ~exist( 'mask', 'var' )
    if D == 1
        mask = ones(Dim,1);
    else
        mask = ones(Dim);
    end
end

boundary = bndry_voxels( logical(mask), 'full' );

% Set the default resadd value
if ~exist('resadd', 'var')
   resadd = 1; 
end

% Set the default do_spm value
if ~exist('do_spm', 'var')
   do_spm = 1; 
end

% Calculate the high resolution mask
mask_hr = mask_highres(mask, resadd);

% xvals = 1:Dim(1);

if D == 1
    Dim = [Dim, 1];
end

if ~exist('lkc_est_version', 'var')
    lkc_est_version = 'conv';
elseif ~(strcmp(lkc_est_version, 'conv') || strcmp(lkc_est_version, 'hpe'))
    error('lkc_est_version must be conv or hpe');
end

%%  Main Function Loop
%--------------------------------------------------------------------------

% Initialize the counters to check the coverage
nabovethresh = 0;
nabovethresh_lat = 0;
nabovethresh_finelat = 0;
nabovethresh_spm = 0;
nabovethresh_lat_spm = 0;
nabovethresh_finelat_spm = 0;

% Main
for b = 1:niters
    %Display b if mod(b,100) = 0
    modul(b,100)
    
    % Calculate data
    lat_data = spfn(sample_size);
    
    % Calculate the convolution t-field and the smooth data. Note this will
    % change later because you'll get this from LKC_est as an output there
    % too!
    [ tfield_finelat, xvals_vecs, smooth_data ] = convfield_t( lat_data.*mask, Kernel, resadd );
    
    % Define the locations of the original voxels (should remove this from
    % the for loop for speed
    voxel_locs = {(ceil(resadd/2) + 1):(resadd+1):length(xvals_vecs{1})};
    orig_lattice_locs = repmat(voxel_locs, 1, D);
    
    % Find the maximum on the lattice
    tfield_lat = tfield_finelat(orig_lattice_locs{:});
    max_tfield_lat = max(tfield_lat(:).*zero2nan(mask(:)));
    
    % Find the maximum on the fine lattice
    max_tfield_finelat = max(tfield_finelat(:).*zero2nan(mask_hr(:)));
    
    % Calculate the LKCs
    if strcmp(lkc_est_version, 'conv')
        LKCs = LKC_conv_est( lat_data, mask, Kernel, resadd );
    elseif strcmp(lkc_est_version, 'hpe')
        % Need to check this is still the right format for the HPE estimate!
        LKCs  = LKC_HP_est( smooth_data, mask_hr, 1 );
    end
    
    % Convert the LKCs to resels
    resel_vec = LKC2resel(LKCs);
    
    % Calculate the RFT threshold using SPM
    threshold = spm_uc_RF_mod(0.05,[1,sample_size-1],'T',resel_vec,1);
    
    % Find the top 3 values on the lattice for initialization
%     [ peak_est_locs, ~, peakvals ] = lmindices(tfield_lat, 3, mask);
    
    % Need to work out the fine lattice evaluation!
    [ peak_est_locs, ~, peakvals ] = lmindices(tfield_finelat, 3, mask_hr);
    
    % Ensure that peak_est_locs is D by npeaks
    if D == 1
        peak_est_locs = peak_est_locs';
    end
    
    % Evaluate the peak values at the xval coordinates
    peak_est_locs = xvaleval(peak_est_locs, xvals_vecs);
    
    % Find the peaks that are within gap of the threshold
    %Need to incorporate the derivative to choose this gap correctly here
    if D == 3
        gap = 5;
    else
        gap = 0.5;
    end
    peaksest2use = peak_est_locs(:, peakvals > (threshold - gap));
    
    % Convert to a cell array in 1D for input to findconvpeaks
    if D == 1
        peaksest2use = num2cell(peaksest2use);
    end
    
    % Find the peaks of the convolution t-field
    if ~isempty(peaksest2use)
        [~, max_tfield_at_lms] = findconvpeaks(lat_data, Kernel, peaksest2use, 'T', mask);
    else
        max_tfield_at_lms = max_tfield_finelat;
    end
    max_tfield_at_lms = max(max_tfield_at_lms);
    
    % This line is included because in 1D convfield and applyconvfield.
    % These are small but need to be resolved. comment this out and
    % uncomment lines 159-161 to see this problem. We need to resolve these
    % differences
    if D == 1
        max_tfield_at_lms = max(max_tfield_at_lms, max_tfield_finelat); 
    end
    if (max_tfield_at_lms + 10^(-10))< max_tfield_finelat
        max_tfield_at_lms
        max_tfield_finelat
        warning('Max lm higher than max lat');
    end

    if max_tfield_at_lms > threshold
        nabovethresh = nabovethresh + 1;
    end
    if  max_tfield_lat > threshold
        nabovethresh_lat = nabovethresh_lat + 1;
    end
    if  max_tfield_finelat > threshold
        nabovethresh_finelat = nabovethresh_finelat + 1;
    end
    
%     max_tfield_lat
%     max_tfield_finelat
%     max_tfield_at_lms
%     pause
%     threshold
%     pause
    
    if do_spm
        orig_lattice_locs{D+1} = ':';
        smooth_data_lat = smooth_data(orig_lattice_locs{:});
        
        % Calculate the resels under the assumption of stationarity
        if D == 3
            lat_FWHM_est = est_smooth(smooth_data_lat, mask); %At the moment this is a biased estimate of course, could use a better convolution estimate of the FWHM and see how that did might be okay for a convolution field?? Would be intersting to see how well it did!
            spm_resel_vec = spm_resels_vol(double(mask),lat_FWHM_est);
        elseif (D == 2 && isequal( mask, ones(Dim))) || (D == 1 && isequal( mask, ones([Dim, 1])))
            lat_FWHM_est = est_smooth(smooth_data_lat, mask); %At the moment this is a biased estimate of course, could use a better convolution estimate of the FWHM and see how that did might be okay for a convolution field?? Would be intersting to see how well it did!
            spm_resel_vec = spm_resels(lat_FWHM_est,Dim, 'B');
        end
        
        threshold_spm = spm_uc_RF_mod(0.05,[1,sample_size-1],'T',spm_resel_vec,1);
        
        % Since the threshold 
        if isempty(peaksest2use) && (peakvals(1) > (threshold_spm - gap))
            [~, max_tfield_at_lms] = findconvpeaks(lat_data, Kernel, peak_est_locs(:,1), 'T', mask, boundary);
        else
            max_tfield_at_lms = max_tfield_finelat;
        end
        if max_tfield_at_lms > threshold_spm
            nabovethresh_spm = nabovethresh_spm + 1;
        end
        if  max_tfield_lat > threshold_spm
            nabovethresh_lat_spm = nabovethresh_lat_spm + 1;
        end
        if  max_tfield_finelat > threshold_spm
            nabovethresh_finelat_spm = nabovethresh_finelat_spm + 1;
        end
    end
end

coverage.conv = nabovethresh/niters;
coverage.lat =  nabovethresh_lat/niters;
coverage.finelat =  nabovethresh_finelat/niters;

if do_spm
    coverage.convspm = nabovethresh_spm/niters;
    coverage.latspm =  nabovethresh_lat_spm/niters;
    coverage.finelatspm =  nabovethresh_finelat_spm/niters;
end

end

