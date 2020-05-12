function [ fwhm_est, Lambda_est, sigma_est ] = est_smooth( data, mask )
% EST_SMOOTH estimates the smoothness of a process. NEED TO DO THE
% CROSS TERMS!! Would be good to study exactly how biased this is.
%--------------------------------------------------------------------------
% ARGUMENTS
% data      Dim by nsubj, data which has nan where there is missing data.
% mask      A mask of the data which is made of 1s and 0s. 1s for where
%           there is data and nans for where there is no data. Default is
%           taken to be the mask with 1s everywhere.
%--------------------------------------------------------------------------
% OUTPUT
% fwhm_est      An estimate of the fwhm in each of the directions.
% Lambda_est    An estimate of the covariance matrix of the partial
%               derivatives.
% sigma_est     An estimate of the smoothness in terms of sigma.
%--------------------------------------------------------------------------
% EXAMPLES
% Dim = 160;
% nsubj = 20;
% noise = noisegen(Dim, nsubj, 6);
% est_smooth(noise)
%
% Dim = [250,250];
% nsubj = 1;
% noise = noisegen(Dim, nsubj, 20);
% est_smooth(noise)
% 
% Dim = [250,250];
% nsubj = 20;
% noise = noisegen(Dim, nsubj, 4);
% est_smooth(noise)
% 
% Dim = [250,250];
% nsubj = 100;
% noise = noisegen(Dim, nsubj, 6);
% est_smooth(noise)
% 
% Dim = [91,109,91];
% nsubj = 100;
% noise = noisegen(Dim, nsubj, 6);
% est_smooth(noise) %Gets around 6.14
% 
% % 3D example with a mask
% Dim = [91,109,91];
% nsubj = 20;
% noise = noisegen(Dim, nsubj, 3);
% mask = imgload('MNImask');
% est_smooth(noise, mask) %Gets around 3.26
%--------------------------------------------------------------------------
% AUTHORS: Samuel Davenport and Fabian Telschow
%--------------------------------------------------------------------------

%% Setup
size_of_data = size(data);
if length(size_of_data) < 5
    nDim = length(size_of_data) - 1;
else
    error('Est_smooth only works in 1,2 and 3 dimensions')
end
% nVox = prod(size_of_data(1:nDim));
nsubj = size_of_data(end);


%% Masking 
%Set the default mask to be all ones.
if nargin < 2
    if nDim == 1
        mask = ones(size_of_data(1), 1);
    else
        mask = ones(size_of_data(1:nDim));
    end
end

%Make the zero-entries of the mask nan.
mask = zero2nan(mask);

% Index to remove cases for different dimensions
index = repmat( {':'}, 1, nDim );
%% Standardize and Mask
mean_over_subjects = mean(data, length(size_of_data) );
nVox = sum(~isnan(mask(:)));

%Subtract the mean and multiply by the mask.
for I = 1:nsubj
    data(index{:}, I) = (data(index{:}, I) - mean_over_subjects).*mask;
end

var_est = sum(data(~isnan(data)).^2)/(nVox*(nsubj - 1));

data = data/sqrt(var_est);

%% Estimate Lambda Matrix and FWHMs
Lambda_est = zeros(nDim);
fwhm_est = zeros(1,nDim);

Xderivmate = diff(data,1,1);
tmp = ~isnan(Xderivmate(index{:}, 1));
denom = sum(tmp(:))*(nsubj-1); %The first half of this
% is the number of voxels in a given subject where Xderivmate is not nan.
%Since we have nsubj subjects we have to multiply this by nsubj - 1 to get
%the total number of voxels. Note its nsubj - 1 rather than nsubj since by
%subtracting the mean we have induced dependence across subjects.
%This is based off of Worsley's 1992 brain paper. Need to check how the
%derivation works to make things unbiased.

Lambda_est(1,1) = sum(Xderivmate(~isnan(Xderivmate)).^2)/denom;
fwhm_est(1) = sqrt(4*log(2)/Lambda_est(1,1));

if nDim > 1
    Yderivmate = diff(data,1,2);
    tmp = ~isnan(Yderivmate(index{:}, 1));
    denom = sum(tmp(:))*(nsubj-1); %The number of non-nans
    Lambda_est(2,2) = sum(Yderivmate(~isnan(Yderivmate)).^2)/denom;
    fwhm_est(2) = sqrt(4*log(2)/Lambda_est(2,2));
end
if nDim > 2
    Zderivmate = diff(data,1,3);
    tmp = ~isnan(Zderivmate(index{:}, 1));
    denom = sum(tmp(:))*(nsubj-1);  
    Lambda_est(3,3) = sum(Zderivmate(~isnan(Zderivmate)).^2)/denom;
    fwhm_est(3) = sqrt(4*log(2)/Lambda_est(3,3));
end

sigma_est = fwhm_est/(sqrt(8*log(2)));

end

%     nZvox = sum(~isnan(Zderivmate(:)));
%     Lambda_est(3,3) = sum(Zderivmate(:).^2)/((nVox-1)*(nsubj-1));
%     Lambda_est(3,3) = sum(Zderivmate(~isnan(Zderivmate)).^2)/((nVox-1)*(nsubj-1));
