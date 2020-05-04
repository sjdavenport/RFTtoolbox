function coverage = record_coverage( spfn, sample_size, FWHM, mask, niters, usehpe )
% RECORD_COVERAGE( data, FWHM, mask, B, sample_size ) estimates the coverage
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
% niters        the number of resamples of the data to do
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
    niters = 1000;
end

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

if D == 3
    usehpe = 1;
else
    usehpe = 0;
    Gker_param = FWHM2sigma(FWHM);
    resadd = 100;
end

nabovethresh = 0;
nabovethresh_lat = 0;
nabovethresh_spm = 0;
nabovethresh_lat_spm = 0;
for b = 1:niters
    b
    boot_lat_data = spfn(sample_size);
    
    [ boot_smoothtfield_lat, ~, smoothed_boot_data ] = smoothtstat( boot_lat_data, FWHM );
    
%     tcf = @(x) tcfield( x, boot_lat_data, FWHM, -1, xvals, mask );
    tcf = @(x) tcfield( x, boot_lat_data, FWHM, -1, xvals, mask );

    
    if usehpe
        HPE  = LKCestim_HPE( smoothed_boot_data, D, mask, 1 );
        L = HPE.hatn;
    else
        L = LKC_GaussConv( boot_lat_data, Gker_param, D, resadd );
    end
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
    elseif D == 1 && isequal( mask, ones([1, Dim])) 
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
%     if max_tfield_at_lms < 0
%         disp('wait')
%     end
    
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

coverage.convhpe = nabovethresh/niters;
coverage.lathpe =  nabovethresh_lat/niters;

if do_spm
    coverage.convspm = nabovethresh_spm/niters;
    coverage.latspm =  nabovethresh_lat_spm/niters;
end

end
