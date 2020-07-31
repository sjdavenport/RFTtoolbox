function coverage = record_coverage_dep( spfn, sample_size, Kernel, resadd, niters, lkc_est_version, do_spm )
% RECORD_COVERAGE( data, FWHM, mask, B, sample_size ) estimates the coverage
% provided by a variety of RFT implementations including non-stationary and
% stationary convolution and lattice versions.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  spfn         a function handle
%  sample_size   the size of each sample to be sampled from the data
%  FWHM          the applied FWHM of the Gaussian Kernel in each direction
%          (we smooth with an istropic Kernel as is commonly done in practice)
% Optional
%  resadd       a non-negative integer giving the resolution increase.
%               Default is 1.
%  niters        the number of resamples of the data to do
%  lkc_est_version      either 'conv' or 'hpe'. Default is 'conv'
%  do_spm       additionally calculate the lkcs using SPM (i.e. under
%               stationarity
%--------------------------------------------------------------------------
% OUTPUT
%  coverage    a structural array with entries:
%     .conv    the average coverage obtained over the number of niters
%               using the convolution approach
%     .lat     coverage obtained via evaluation on a lattice (this is
%              conservative)
%     .finelat   the coverage obtained with evaluation on the fine lattice
%              given by resadd, this is conservative but not as bad as 
%             coverage.lat. On small domains or for very smooth fields this
%             will be similar to coverage.conv, especially for high resadd.
%--------------------------------------------------------------------------
% EXAMPLES
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
single_sample_field = spfn(1);
if strcmp(class(single_sample_field), 'Field')
    direct_field = 1;
else
    single_sample_field = single_sample_field.lat_data;
    direct_field = 0;
end
Dim = single_sample_field.fieldsize;
D = single_sample_field.D;
mask = single_sample_field.mask;

% Calculate the Euler characteristic
L0 = EulerChar(mask, 0.5, D);

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'niters', 'var' )
    % default option of niters
    niters = 1000;
end

% Set the default resadd value
if ~exist('resadd', 'var')
   resadd = 1; 
end

% Set the default do_spm value
if ~exist('do_spm', 'var')
   do_spm = 1; 
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
    if direct_field == 1
        lat_data = spfn(sample_size);
    else
        lat_data = spfn(sample_size).lat_data;
    end
    
    % Calculate the fine t convolution field and the individual convolution
    % fields
    [ tcfield, cfields ] = convfield_t_Field( lat_data, Kernel, resadd );
    
    % Calculate the derivatives of the convolution fields
    dcfields = convfield_Field( lat_data, Kernel, 1, resadd, 1 );

    % Define the locations of the original voxels
    orig_lattice_locs = cell(1,D);
    for d = 1:D
        orig_lattice_locs{d} = (ceil(resadd/2) + 1):(resadd+1):length(cfields.xvals{d});
    end
    
    % Find the maximum on the lattice
    tfield_lat = tcfield.field(orig_lattice_locs{:});
    max_tfield_lat = max(tfield_lat(:).*zero2nan(mask(:)));
    
    % Find the maximum on the fine lattice
    max_tfield_finelat = max(tcfield.field(:).*zero2nan(cfields.mask(:)));
    
    % Calculate the LKCs
    if strcmp(lkc_est_version, 'conv')
        % Obtain the LKCs using the convolution estimate
        L = LKC_voxmfd_est( cfields, dcfields ); %L0 in LKC_voxmfd_est calculated wrong: using the highres mask
%         LKCs = LKC_conv_est( lat_data, mask, Kernel, resadd );
    elseif strcmp(lkc_est_version, 'hpe')
        % Obtain the LKCs using the HPE
        LKCs  = LKC_HP_est( cfields, 1, 1 );
        L = LKCs.hatL; L0 = LKCs.L0;
    end
    
    % Convert the LKCs to resels
    resel_vec = LKC2resel(L, L0);
    
    % Calculate the RFT threshold using SPM
    threshold = spm_uc_RF_mod(0.05,[1,sample_size-1],'T',resel_vec,1);
    
    % Initialize peak_est_locs and peakvals
    [ peak_est_locs, ~, peakvals ] = lmindices(tcfield.field, 3, cfields.mask);
    
    % Ensure that peak_est_locs is D by npeaks
    if D == 1
        peak_est_locs = peak_est_locs';
    end
    
    % Evaluate the peak values at the xval coordinates
    peak_est_locs = xvaleval(peak_est_locs, cfields.xvals);
    
    % Find the peaks that are within gap of the threshold
    % Need to incorporate the derivative to choose this gap correctly here
    if D == 3
        gap = 5;
    else
        gap = 1;
    end
    peaksest2use = peak_est_locs(:, peakvals > (threshold - gap));
    
    % Convert to a cell array in 1D for input to findconvpeaks
    if D == 1
        peaksest2use = num2cell(peaksest2use);
    end
    
    % Find the peaks of the convolution t-field
    if ~isempty(peaksest2use)
        [~, max_tfield_at_lms] = findconvpeaks(lat_data.field, Kernel, peaksest2use, 'T', mask);
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
        smooth_data_lat = cfields.field(orig_lattice_locs{:});
        
        % Calculate the resels under the assumption of stationarity
        if D == 3
            lat_FWHM_est = est_smooth(smooth_data_lat, mask); %At the moment this is a biased estimate of course, could use a better convolution estimate of the FWHM and see how that did might be okay for a convolution field?? Would be intersting to see how well it did!
            spm_resel_vec = spm_resels_vol(double(mask),lat_FWHM_est);
        elseif (D == 2 && isequal( mask, ones(Dim))) || (D == 1 && isequal( mask, ones([Dim, 1])))
            warning('This only works for a box!')
            lat_FWHM_est = est_smooth(smooth_data_lat, mask); %At the moment this is a biased estimate of course, could use a better convolution estimate of the FWHM and see how that did might be okay for a convolution field?? Would be intersting to see how well it did!
            spm_resel_vec = spm_resels(lat_FWHM_est,Dim, 'B');
        end
        
        threshold_spm = spm_uc_RF_mod(0.05,[1,sample_size-1],'T',spm_resel_vec,1);
        
        % Since the threshold 
        if isempty(peaksest2use) && (peakvals(1) > (threshold_spm - gap))
            [~, max_tfield_at_lms] = findconvpeaks(lat_data.field, Kernel, peak_est_locs(:,1), 'T', mask);
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

