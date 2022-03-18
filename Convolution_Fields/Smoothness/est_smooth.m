function [ fwhm_est_forman, fwhm_est_kiebel, Lambda_est, sigma_est, ...
    fwhm_est_forman_unscaled, fwhm_est_kiebel_unscaled, Lambda_est_unscaled]...
    = est_smooth( data, combine_var, mask, lat_spacing, df )
% EST_SMOOTH estimates the smoothness of a process. NEED TO DO THE
% CROSS TERMS!! Would be good to study exactly how biased this is.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  data      Dim by nsubj, data which has nan where there is missing data.
% Optional
%  mask      A mask of the data which is made of 1s and 0s. 1s for where
%           there is data and nans for where there is no data. Default is
%           taken to be the mask with 1s everywhere.
%  combine_var  0/1 whether or not to combine the pool across all subjects
%              to estimate the variance. Only valid if the data has the
%              same variance at each voxel (so not for fMRI data!). If
%              combine_var is -1, then no variance scaling is done at all.
%  lat_spacing   the distance between adjacent voxels. Default is 1. Taking
%                smaller values allows for fine lattice input.
%--------------------------------------------------------------------------
% OUTPUT
% fwhm_est      An estimate of the fwhm in each of the directions.
% Lambda_est    An estimate of the covariance matrix of the partial
%               derivatives.
% sigma_est     An estimate of the smoothness in terms of sigma.
%--------------------------------------------------------------------------
% EXAMPLES
% Dim = 160; nsubj = 20;
% noise = noisegen(Dim, nsubj, 2);
% % Without variance pooling
% [ fwhm_est_forman, fwhm_est_kiebel] = est_smooth(noise);
% % With variance pooling
% [ fwhm_est_forman_pv, fwhm_est_kiebel_pv] = est_smooth(noise, ones(Dim, 1), 1);
% fwhm_est_forman, fwhm_est_kiebel, fwhm_est_forman_pv, fwhm_est_kiebel_pv
%
% Dim = [50,50];
% nsubj = 1;
% noise = noisegen(Dim, nsubj, 20);
% est_smooth(noise)
%
% Dim = [250,250];
% nsubj = 100;
% noise = noisegen(Dim, nsubj, 20);
% [ fwhm_est_forman, fwhm_est_kiebel]  = est_smooth(noise)
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

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist('combine_var', 'var')
    combine_var = 0;
end

if ~exist('lat_spacing', 'var')
    lat_spacing = 1;
end

%% Setup
size_of_data = size(data);
if length(size_of_data) < 5
    D = length(size_of_data) - 1;
else
    error('Est_smooth only works in 1,2 and 3 dimensions')
end
% nVox = prod(size_of_data(1:nDim));
nsubj = size_of_data(end);

%% Masking
%Set the default mask to be all ones.
if ~exist('mask', 'var')
    if D == 1
        mask = ones(size_of_data(1), 1);
    else
        mask = ones(size_of_data(1:D));
    end
end

%Make the zero-entries of the mask nan.
mask = zero2nan(mask);

% Index to remove cases for different dimensions
index = repmat( {':'}, 1, D );

% Compute the number of voxels within the mask
nVox  = sum( ~isnan( mask(:) ) );

%% Standardize and Mask
if ~exist('df', 'var')
    mean_over_subjects = mean( data, length( size_of_data ) );
    
    % Subtract the mean and multiply by the mask.
    for I = 1:nsubj
        data(index{:}, I) = (data(index{:}, I) - mean_over_subjects).*mask;
    end
    df = 1;
else
    % Mask the data
    data = data .* mask;
end

if combine_var > -1 %Outer loop allows for removing variance scaling completely
    if combine_var  % Assumes that the variance is the same at each voxel
        % and pools over voxels to estimate it
        var_est = sum(data(~isnan(data)).^2)/(nVox*(nsubj - 1));
        data = data/sqrt(var_est);
    else % Assumes that the variance is different at each voxel
        std_dev = std(data, 0, D + 1);
        data = data./std_dev;
    end
end

% Allow for estimation where the lattice is finer than 1 point
data = data/lat_spacing;

%% Estimate Lambda Matrix and FWHMs
Lambda_est = zeros(D);
fwhm_est_kiebel = zeros(1,D);

Xderivmate = diff(data,1,1);
Xdirectionmask = ~isnan(Xderivmate(index{:}, 1));

denom = obtain_denom(Xdirectionmask, nsubj, combine_var);
Lambda_est(1,1) = sum(Xderivmate(~isnan(Xderivmate)).^2)/denom;
fwhm_est_kiebel(1) = sqrt(4*log(2)/Lambda_est(1,1));

if D > 1
    Yderivmate = diff(data,1,2);
    Ydirectionmask = ~isnan(Yderivmate(index{:}, 1));
    denom = obtain_denom(Ydirectionmask, nsubj, combine_var);
    Lambda_est(2,2) = sum(Yderivmate(~isnan(Yderivmate)).^2)/denom;
    fwhm_est_kiebel(2) = sqrt(4*log(2)/Lambda_est(2,2));
end
if D > 2
    Zderivmate = diff(data,1,3);
    Zdirectionmask = ~isnan(Zderivmate(index{:}, 1));
    denom = obtain_denom(Zdirectionmask, nsubj, combine_var);
    Lambda_est(3,3) = sum(Zderivmate(~isnan(Zderivmate)).^2)/denom;
    fwhm_est_kiebel(3) = sqrt(4*log(2)/Lambda_est(3,3));
end

sigma_est = fwhm_est_kiebel/(sqrt(8*log(2)));

sigma_est_forman = abs(sqrt(-1./4./log(1-diag(Lambda_est)/2))); % This can be negative for low FWHM
fwhm_est_forman = sigma2FWHM(sigma_est_forman);

if ~combine_var
    Lambda_est_unscaled = Lambda_est*(nsubj-2)/(nsubj-3);
    fwhm_est_kiebel_unscaled = sqrt(4*log(2)./diag(Lambda_est_unscaled));
    sigma_est_forman_unscaled = abs(sqrt(-1./4./log(1-diag(Lambda_est_unscaled)/2))); % This can be negative for low FWHM
    fwhm_est_forman_unscaled = sigma2FWHM(sigma_est_forman_unscaled);
else
    fwhm_est_kiebel_unscaled = NaN;
    fwhm_est_forman_unscaled = NaN;
    Lambda_est_unscaled = NaN;
end

end

function denom = obtain_denom(directional_mask, nsubj, combine_var)
if combine_var
    denom = sum(directional_mask(:))*(nsubj-1);
else
    % Include the scaling factor for unbiased estimation
    denom = sum(directional_mask(:))*(nsubj-1)*(nsubj-2)/(nsubj-3);
    %     denom = sum(tmp(:))*(nsubj-1);
end
% denom = sum(tmp(:))*(nsubj-df); %The first half of this
% is the number of voxels in a given subject where Xderivmate is not nan.
%Since we have nsubj subjects we have to multiply this by nsubj - 1 to get
%the total number of voxels. Note its nsubj - 1 rather than nsubj since by
%subtracting the mean we have induced dependence across subjects.
%This is based off of Worsley's 1992 brain paper. Need to check how the
%derivation works to make things unbiased.
end