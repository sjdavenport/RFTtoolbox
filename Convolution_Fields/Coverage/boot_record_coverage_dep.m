function coverage = record_coverage_boot( lat_data, sample_size, FWHM, mask, B )
% RECORD_COVERAGE_BOOT( data, FWHM, mask, B, sample_size ) estimates the coverage
% provided by a variety of RFT implementations including non-stationary and
% stationary convolution and lattice versions.
%--------------------------------------------------------------------------
% ARGUMENTS
% data          a Dim by sample_size array of data. D = length(Dim) must be <= 3
% sample_size   the size of each sample to be sampled from the data
% FWHM          the applied FWHM of the Gaussian Kernel in each direction
%          (we smooth with an istropic Kernel as is commonly done in practice)
% mask          a 0/1 array of size Dim which provides a mask of the data
%               the default to use is no mask i.e. 1
% B             the number of resamples (i.e. bootstraps) of the data to do
%--------------------------------------------------------------------------
% OUTPUT
% coverage is a structure giving the coverage for each setting
% coverage.convhpe   HPE convolution coverage
% coverage.lathpe    HPE lattice coverage
% coverage.convspm   SPM (using FWHM estimation) convolution coverage
% coverage.lathpe    SPM (using FWHM estimation) lattice coverage
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if nargin < 5
    B = 1000;
end

if ischar(lat_data)
    % this is if you need to continuously load in files because they are too
    % large to fit all of them in memory at once.
    if strcmp(lat_data, 'RS_2Block')
        global featrunloc
        imgdir = [featrunloc, 'RS_2Block', '_warped/'];
    else
        error('Not set up for this data input')
    end
    totalnoofsubs = length(filesindir(imgdir));
    Sdata = size(lat_data);
    totalnoofsubs = Sdata(end);
    Dim = Sdata(1:end-1);
    D = length(Dim);
else
    Sdata = size(lat_data);
    totalnoofsubs = Sdata(end);
    Dim = Sdata(1:end-1);
    D = length(Dim);
    [~,~,smooth_data] = smoothtstat(lat_data,FWHM);
end

if nargin < 4
    if D == 1
        mask = ones(1,Dim);
    else
        mask = ones(Dim);
    end
end

xvals = 1:Dim(1);

if D == 1
    Dim = [1,Dim];
end

nabovethresh = 0;
nabovethresh_lat = 0;
nabovethresh_spm = 0;
nabovethresh_lat_spm = 0;
for b = 1:B
    b
    sample_index = randsample(totalnoofsubs,sample_size);
    
    if isstr(lat_data)
        boot_lat_data = zeros([91,109,91,sample_size]);
        
        for K = 1:length(sample_index)
%             boot_lat_data(:,:,:,K) = readRSD( K );
        end
        
        [ boot_smoothtfield_lat, ~, smoothed_boot_data ] = smoothtstat( boot_lat_data, FWHM );
    else
        if D == 1
            boot_lat_data = lat_data(:, sample_index);
            smoothed_boot_data = smooth_data(:, sample_index);
        elseif D == 2
            boot_lat_data = lat_data(:,:, sample_index);
            smoothed_boot_data = smooth_data(:,:, sample_index);
        elseif D == 3
            boot_lat_data = lat_data(:,:,:, sample_index);
            smoothed_boot_data = smooth_data(:,:,:, sample_index);
        end
        boot_smoothtfield_lat = mvtstat(smoothed_boot_data, Dim);
    end
    
    tcf = @(x) tcfield( x, boot_lat_data, FWHM, -1, xvals, mask );
    
    HPE  = LKCestim_HPE( smoothed_boot_data, D, mask, 1 );
    L = HPE.hatn;
    %     resel_vec = [1, L(1)/sqrt(4*log(2)),  L(2)/(4*log(2)), L(3)/(4*log(2))^(3/2)];
    resel_vec = LKC2resel(L);
    threshold = spm_uc_RF_mod(0.05,[1,sample_size-1],'T',resel_vec,1);
    
    if D == 3
        lat_FWHM_est = est_smooth(smoothed_boot_data, mask); %At the moment this is a biased estimate of course, could use a better convolution estimate of the FWHM and see how that did might be okay for a convolution field?? Would be intersting to see how well it did!
        spm_resel_vec = spm_resels_vol(mask,lat_FWHM_est);
        do_spm = 1;
    elseif D == 2 && isequal( mask, ones(Dim) )
        lat_FWHM_est = est_smooth(smoothed_boot_data, mask); %At the moment this is a biased estimate of course, could use a better convolution estimate of the FWHM and see how that did might be okay for a convolution field?? Would be intersting to see how well it did!
        spm_resel_vec = spm_resels(lat_FWHM_est,Dim, 'B');
        do_spm = 1;
    else
        do_spm = 0;
    end
    
    peak_est_locs = lmindices(boot_smoothtfield_lat, 3, mask);
    
    top_lmlocs = findconvpeaks_t(boot_lat_data, FWHM, peak_est_locs, mask);
    
    max_tfield_lat = max(boot_smoothtfield_lat(:).*zero2nan(mask(:)));
    max_tfield_at_lms = max(tcf(top_lmlocs));
    max_tfield_at_lms = max(max_tfield_lat, max_tfield_at_lms);
    
    if max_tfield_at_lms > threshold
        nabovethresh = nabovethresh + 1;
    end
    if  max_tfield_lat > threshold
        nabovethresh_lat = nabovethresh_lat + 1;
    end
    
    if do_spm
        threshold_spm = spm_uc_RF_mod(0.05,[1,sample_size-1],'T',spm_resel_vec,1);
        
        if max_tfield_at_lms > threshold_spm
            nabovethresh_spm = nabovethresh + 1;
        end
        if  max_tfield_lat > threshold_spm
            nabovethresh_lat_spm = nabovethresh_lat + 1;
        end
    end
end

coverage.convhpe = nabovethresh/B;
coverage.lathpe =  nabovethresh_lat/B;

if do_spm
    coverage.convspm = nabovethresh_spm/B;
    coverage.latspm =  nabovethresh_lat_spm/B;
end

end
