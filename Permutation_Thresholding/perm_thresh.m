function [im_perm, threshold, vec_of_maxima] = perm_thresh( data, stat, FWHM, mask, demean, niters, include_original, alpha )
% PERM_THRESH( data, stat, niters, include_original, subject_mask, alpha ) 
% implements the permutation test voxelwise on a dataset in order to 
% estimate a one-sample threshold with which to perform multiple 
% testing. This makes the assumption that the errors are symmetric.
%--------------------------------------------------------------------------
% ARGUMENTS
% data      a Dim by nsubj array of the data.
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
% data = normrnd(0,1,20,100);
% [~, threshold] = perm_thresh(data, 'T');
%
% lat_data = wfield([5,5], 10);
% tic; [~,threshold] = perm_thresh_new(lat_data, 'T'); toc
% tic; [~,threshold] = perm_thresh(lat_data.field, 'T');toc
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
if ~exist('subject_mask', 'var')
    mask = NaN;
end
if ~exist('demean', 'var')
    demean = 0;
end
if ~exist('alpha', 'var')
    alpha = 0.05;
end

% Obtain the size of the data input
sD = size(data);

% Calculate the number of subjects
nsubj = sD(end);

% Calculate the number of dimensions
D = length(sD) - 1;
    
if ~isnan(FWHM)
    data = fconv(data, FWHM, D);
end

% Obtain the original tstat
if strcmp(stat, 'Z')
    stat_image = mean(data, D+1);
elseif strcmp(stat, 'T')
    stat_image = mvtstat(data);
end
% 
if demean
   data = data - mean(data, D+1); 
end

if include_original
%     niters = niters - 1;
    start = 2;
    
    if strcmp(stat, 'Z')
        combined_orig = mean(data, (D+1));
    elseif strcmp(stat, 'T')
        combined_orig = mvtstat(data);
    end
    
    vec_of_maxima = zeros(1, niters );
    if isnan(mask)
        vec_of_maxima(1) = max(combined_orig(:));
    else
        vec_of_maxima(1) = combined_orig(lmindices(combined_orig, 1, mask));
    end
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
    modul(iter, 100);
    random_berns_for_iter = random_berns(:, iter);
    random_sample = find(random_berns_for_iter < 0);
    data2 = data;
    data2(indexD{:}, random_sample) = -data(indexD{:}, random_sample); %This doesn't subtract the mean though it probably should!
    
    if strcmp(stat, 'Z')
        combined = mean(data2, (D+1));
    elseif strcmp(stat, 'T')
        combined = mvtstat( data2 );
    end
    
    if isnan(mask)
        vec_of_maxima(iter) = max(combined(:));
    else
%         vec_of_maxima(iter) = max(combined(:).*mask(:));
        vec_of_maxima(iter) = combined(lmindices(combined, 1, mask));
    end
end

threshold = prctile(vec_of_maxima, 100*(1-alpha) );
im_perm = stat_image > threshold;

end