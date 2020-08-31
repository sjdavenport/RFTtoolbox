 function coverage = record_coverage( spfn, sample_size, params, niters, npeaks, version )
% RECORD_COVERAGE( data, FWHM, mask, B, sample_size ) estimates the coverage
% provided by a variety of RFT implementations including non-stationary and
% stationary convolution and lattice versions.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  spfn         a function handle which given an integer number of subjects
%               nsubj generates nsubj random fields and saves them in a
%               data type of class field
%  sample_size  the size of each sample to be sampled from the data
%  params       an object of class ConvFieldParams
% Optional
%  niters        the number of resamples of the data to do
%  npeaks       the number of peaks of the lattice data around which to
%               search for the local maxima
%  version      the version of the LKC estimation to use
%  subsets
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

% Obtain a sample of the field
single_sample_field = spfn(1);

% Obtain properties of the data
D = single_sample_field.D;

% mask = single_sample_field.mask;
% Calculate the Euler characteristic
% L0 = EulerChar(mask, 0.5, D);

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'niters', 'var' )
    % default option of niters
    niters = 1000;
end

% Set the default do_spm value
if ~exist('version', 'var')
    if D < 3
        version = true(ones(1,D));
    elseif D == 3
        version = [true, true, false];
    end
end

if ~exist('npeaks', 'var')
    npeaks = 3;
end

% Determine whether to use fixed random subsets in advance instead of a
% given number of subjects
if iscell(sample_size)
    %Obtain the subsets
    subsets = sample_size;
    
    % Determine the sample size to use
    sample_size = length(subsets{1});
    
    % Calculate the number of subsets stored
    sample_niters = length(subsets);
    
    % Ensure that the number of iterations is the same as the number of
    % samples provided
    if sample_niters ~= niters
        error('The number of subsets must be the same as the number of iterations\n')
    end
    
    % Define an indicator that shows you're using random subsets
    use_subsets = 1;
else
    % Define an indicator that shows you're not using random subsets
    use_subsets = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------

% Initialize the counters to check the coverage
nabovethresh = 0;
nabovethresh_lat = 0;
nabovethresh_finelat = 0;

% Initialize vectors to store the maxima
latmaxima     = zeros( 1, niters );
finelatmaxima = zeros( 1, niters );
convmaxima    = zeros( 1, niters );
thresholds    = zeros( 1, niters );
maxabovethreshold = zeros( 1, niters );

storeLKCs = zeros( D, niters );

allmaxima = zeros( npeaks, niters );

% Initialize the coverage structure
coverage = struct();

% Main
for b = 1:niters
    %Display b if mod(b,100) = 0
    modul(b,10)
    
    % Obtain the data
    if use_subsets
        lat_data = spfn(subsets{b});
    else
        lat_data = spfn(sample_size);
    end
    
    lat_data = Mask(lat_data);
    [ ~, threshold, maximum, L ] = vRFT(lat_data, params, npeaks, version);
    storeLKCs(:,b) = L';
    if any(isnan(L))
        warning('NAN LKC recorded')
    end
    
    if maximum.conv > threshold
        nabovethresh = nabovethresh + 1;
    end
    if  maximum.lat > threshold
        nabovethresh_lat = nabovethresh_lat + 1;
    end
    if  maximum.finelat > threshold
        nabovethresh_finelat = nabovethresh_finelat + 1;
    end
    latmaxima(b) = maximum.lat;
    finelatmaxima(b) = maximum.finelat;
    convmaxima(b) = maximum.conv;
    allmaxima(1:length(maximum.allmaxima),b) = maximum.allmaxima';
    thresholds(b) = threshold;
    maxabovethreshold(b) = sum( maximum.allmaxima > threshold );
    % Error checking loop 
    if maximum.finelat > maximum.conv + 10^(-2)
        a = 1
    end
end

coverage.conv = nabovethresh/niters;
coverage.lat =  nabovethresh_lat/niters;
coverage.finelat =  nabovethresh_finelat/niters;

coverage.finelatmaxima = finelatmaxima;
coverage.latmaxima  = latmaxima;
coverage.convmaxima = convmaxima;
coverage.allmaxima  = allmaxima;
coverage.thresholds = thresholds;

coverage.maxabovethreshold = maxabovethreshold;

coverage.storeLKCs = storeLKCs;

end

% % 
% nabovethresh_spm = 0;
% nabovethresh_lat_spm = 0;
% nabovethresh_finelat_spm = 0;
% % The spm loop needs redoing!
% if do_spm
%     orig_lattice_locs{D+1} = ':';
%     smooth_data_lat = cfields.field(orig_lattice_locs{:});
%     
%     % Calculate the resels under the assumption of stationarity
%     if D == 3
%         lat_FWHM_est = est_smooth(smooth_data_lat, mask); %At the moment this is a biased estimate of course, could use a better convolution estimate of the FWHM and see how that did might be okay for a convolution field?? Would be intersting to see how well it did!
%         spm_resel_vec = spm_resels_vol(double(mask),lat_FWHM_est);
%     elseif (D == 2 && isequal( mask, ones(Dim))) || (D == 1 && isequal( mask, ones([Dim, 1])))
%         warning('This only works for a box!')
%         lat_FWHM_est = est_smooth(smooth_data_lat, mask); %At the moment this is a biased estimate of course, could use a better convolution estimate of the FWHM and see how that did might be okay for a convolution field?? Would be intersting to see how well it did!
%         spm_resel_vec = spm_resels(lat_FWHM_est,Dim, 'B');
%     end
%     
%     threshold_spm = spm_uc_RF_mod(0.05,[1,sample_size-1],'T',spm_resel_vec,1);
%     
%     % Since the threshold
%     if isempty(peak_est_locs) && (peakvals(1) > (threshold_spm - gap))
%         [~, max_tfield_at_lms] = findconvpeaks(lat_data.field, Kernel, peak_est_locs(:,1), 'T', mask);
%     else
%         max_tfield_at_lms = max_tfield_finelat;
%     end
%     if max_tfield_at_lms > threshold_spm
%         nabovethresh_spm = nabovethresh_spm + 1;
%     end
%     if  max_tfield_lat > threshold_spm
%         nabovethresh_lat_spm = nabovethresh_lat_spm + 1;
%     end
%     if  max_tfield_finelat > threshold_spm
%         nabovethresh_finelat_spm = nabovethresh_finelat_spm + 1;
%     end
% end
% if do_spm
%     coverage.convspm = nabovethresh_spm/niters;
%     coverage.latspm =  nabovethresh_lat_spm/niters;
%     coverage.finelatspm =  nabovethresh_finelat_spm/niters;
% end