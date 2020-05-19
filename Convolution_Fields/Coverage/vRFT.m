function [ output_image, threshold, max_finelat ] = vRFT( lat_data, FWHM, spacing, alpha )
% vRFT( lat_data, FWHM, spacing, alpha ) runs voxelwise RFT inference on a set of images to
% detect areas of activation. (no mask implemented yet!)
%--------------------------------------------------------------------------
% ARGUMENTS
% lat_data      a Dim by nsubj array of data
% FWHM          the FWHM with which to smooth
% alpha         the alpha level at which to threshold. Default is 0.05.
%               Recommend alpha <= 0.05 for best performance.
%--------------------------------------------------------------------------
% OUTPUT
% output_image  the (fine lattice) output image
% threshold     the voxelwise RFT threshold
% max_finelat   the maximum on a fine lattice given by spacing
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
sld = size(lat_data);
D = length(size(lat_data)) - 1;
if nargin < 3
    if D == 1 || D == 2
        spacing = 0.25;
    elseif D == 3
        spacing = 0.25;
    end
end
if D == 3
    error('Not working in D = 3 yet')
end
if nargin < 4
    alpha = 0.05;
end
% calculate the convolution field on a finer lattice

tfield_fine = convfield_t(lat_data,FWHM,spacing);

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

max_finelat = max(tfield_fine(:));

% Calculate L and the threshold
resAdd = 2;
L = LKCestim_GaussConv3( lat_data, FWHM, ones(sld(1:end-1)), resAdd);
resel_vec = LKC2resel(L);

% Calculate the threshold
nsubj = size(lat_data,D+1);
threshold = spm_uc_RF(alpha,[1,nsubj-1],'T',resel_vec,1);

output_image = tfield_fine > threshold;
output_image = double(output_image);

end

