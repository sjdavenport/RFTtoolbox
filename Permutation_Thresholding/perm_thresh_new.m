function [im_perm, threshold, vec_of_maxima] = perm_thresh_new( lat_data, stat, FWHM, demean, niters, include_original, alpha )
% PERM_THRESH( data, stat, niters, include_original, subject_mask, alpha ) 
% implements the permutation test voxelwise on a dataset in order to 
% estimate a one-sample threshold with which to perform multiple 
% testing. This makes the assumption that the errors are symmetric.
%--------------------------------------------------------------------------
% ARGUMENTS
% lat_data  lattice data contained within a Field structure
% stat      is the statistic to use. STAT = 'Z', we assume Gaussian and use
%           the mean. STAT = 'T', we compute a one-sample t-stat. Default
%           is STAT = 'Z'.
% niters    the number of permutations to calculate
% include_original 0/1 whether to include the original sample in your
%           niters number of permutations.
% alpha     the alpha level, default is 0.05.
%--------------------------------------------------------------------------
% OUTPUT
% threshold         the threshold calculated from the 100(1-alpha)% percent
%                   quantile of the vector of the maxima
% vec_of_maxima     the 1 by niters vector of the maxima of the random 
%                   lattice fields for each permutation
%--------------------------------------------------------------------------
% EXAMPLES
% lat_data = wfield([10,10], 30);
% [im, threshold] = perm_thresh(lat_data, 'T')
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if ~exist('stat', 'var')
    stat = 'T';
end
if ~exist('niters', 'var')
    niters = 5000;
end
if ~exist('include_original', 'var')
    include_original = 1;
end
if ~exist('FWHM', 'var')
    FWHM = NaN;
end
if ~exist('demean', 'var')
    demean = 0;
end
if ~exist('alpha', 'var')
    alpha = 0.05;
end

% Obtain the size of the data input
sD = size(lat_data);

% Calculate the number of subjects
nsubj = sD(end);

% Calculate the number of dimensions
D = length(sD) - 1;
    
if ~isnan(FWHM)
    resadd = 0;
    params = ConvFieldParams(repmat(FWHM, 1, D), resadd);
    lat_data = convfield(lat_data, params);
end

% Obtain the original tstat
if strcmp(stat, 'Z')
    stat_image = mean(lat_data.field, D+1);
elseif strcmp(stat, 'T')
    stat_image = mvtstat(lat_data.field);
end

if demean
   lat_data.field = lat_data.field - mean(lat_data.field, D+1); 
end

if include_original
    start = 2;
    
    if strcmp(stat, 'Z')
        combined_orig = mean(lat_data.field, (D+1));
    elseif strcmp(stat, 'T')
        combined_orig = mvtstat(lat_data.field);
    end
    
    vec_of_maxima = zeros(1, niters );
    dataonmask = combined_orig(lat_data.mask);
    vec_of_maxima(1) = max(dataonmask(:));
else
    vec_of_maxima = zeros(1, niters);
    start = 1;
end

% Index for multidimensional subsetting
indexD = repmat( {':'}, 1, D );

% Bernoulli random variables for sign flipping
random_berns = 2*(binornd(1,0.5, nsubj, niters )-1/2);

% Main loop
for iter = start:niters
%     modul(iter,1000);
    random_berns_for_iter = random_berns(:, iter);
    random_sample = find(random_berns_for_iter < 0);
    data2 = lat_data.field;
    data2(indexD{:}, random_sample) = -lat_data.field(indexD{:}, random_sample); %This doesn't subtract the mean though it probably should!
    
    if strcmp(stat, 'Z')
        combined = mean(data2, (D+1));
    elseif strcmp(stat, 'T')
        combined = mvtstat( data2 );
    end
    
    dataonmask = combined(lat_data.mask);
    vec_of_maxima(iter) = max(dataonmask(:));
end

threshold = prctile(vec_of_maxima, 100*(1-alpha) );
im_perm = stat_image > threshold;

end