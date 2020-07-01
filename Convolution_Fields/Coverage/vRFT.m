 function [ output_image, threshold, max_finelat, xvals ] = vRFT( lat_data, ...
                                              Kernel, mask, resadd, alpha )
% vRFT( lat_data, Kernel, mask, resadd, alpha ) runs voxelwise RFT 
% inference on a set of images to detect areas of activation using a 
% one-sample t-statistic.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  lat_data      a Dim by nsubj array of data
%  Kernel    either an object of class SepKernel or a numeric.
%            If Kernel is numeric, the convolution fields are generated by 
%            smoothing with an isotropic Gaussian kernel with FWHM = Kernel.
% Optional
%  mask       a logical array where 1s indicate the presence of data.
%             Default is not to mask i.e. to take mask = true(Dim).
%  resadd     the amount of voxels added equidistantly in between the
%             existing voxels. Default is 1.
%  alpha         the alpha level at which to threshold. Default is 0.05.
%                Recommend alpha <= 0.05 for best performance.
%--------------------------------------------------------------------------
% OUTPUT
%  output_image  the (fine lattice) output image
%  threshold     the voxelwise RFT threshold
%  max_finelat   the maximum on a fine lattice given by spacing
%--------------------------------------------------------------------------
% DEVELOPERS TODOS: get rid of the SPM dependence :)
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Obtain the size of the lattice input
s_lat = size(lat_data);

% Obtain the dimensions of the data
Dim = s_lat(1:end-1);

% Obtain the number of dimensions
D = length(Dim);

% Obtain the number of subjects
nsubj = size(lat_data,D+1);

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist('mask', 'var')
    if D == 1
        mask = true(Dim,1);
    else
        mask = true(Dim);
    end
end
if ~exist( 'resadd', 'var' )
    resadd = 1;
end
if ~exist( 'alpha', 'var' )
    alpha = 0.05;
end

%%  Main function
%--------------------------------------------------------------------------
% calculate the convolution t field (uses the default enlarge of
% ceil(resadd/2)) % Need to change so that convY in obtained from
% LKC_conv_est
[tfield_fine, xvals] = convfield_t(lat_data.*mask,Kernel,resadd);
max_finelat = max(tfield_fine(:));

% Calculate L and the threshold
LKCs = LKC_conv_est( lat_data, mask, Kernel, resadd );

% Compute the resels from the LKCs (essentially a scaling factor)
resel_vec = LKC2resel(LKCs);

% Calculate the threshold using SPM (for now) 
threshold = spm_uc_RF(alpha,[1,nsubj-1],'T',resel_vec,1);

% Determine the areas of the image where the t-field exceeds the threshold
output_image = tfield_fine > threshold;

% Convert the output to double (can't remember why I do this!)
output_image = double(output_image);

end

% DEPRECATED
% high_local_maxima = lmindices(tfield_fine, 3);
% Calculate initial estimates of peak location
% if D == 1
%     peak_est_locs = [NaN,setdiff(xvals_fine(high_local_maxima),[1,nvox])];
% end

% tcf = @(x) tcfield( x, lat_data, FWHM );

% if length(peak_est_locs) == 1
%     tfield_at_lms = -Inf; %If the local max occurs at the boundary you don't need to account for it.
% else
%     top_lmlocs = findconvpeaks_t(lat_data, FWHM, peak_est_locs);
%     tfield_at_lms = tcf(top_lmlocs);
% end
% 
% % Calculate the maximum on the lattice and of the convolution field
% max_finelat = max(tfield_fine);
% max_conv = max([tfield_at_lms,max_finelat]); %Included for stability in case the maximum finding didn't work correctly.
