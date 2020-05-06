function [threshold, vec_of_maxima] = perm_thresh( data, stat, niters, include_original, subject_mask, alpha )
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
% data = normrnd(0,1, 100, 5);
% threshold = perm_thresh(data)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport.
if nargin < 2
    stat = 'T';
end
if nargin < 3
    niters = 5000;
end
if nargin < 4
    include_original = 1;
end
if nargin < 5
    subject_mask = NaN;
end
if nargin < 6
    alpha = 0.05;
end

sD = size(data);
nsubj = sD(end);
D = length(sD) - 1;

if include_original
%     niters = niters - 1;
    start = 2;
    
    if strcmp(stat, 'Z')
        combined = mean(data, (D+1));
    elseif strcmp(stat, 'T')
        combined = mvtstat(data);
    end
    
    vec_of_maxima = zeros(1, niters );
    if isnan(subject_mask)
        vec_of_maxima(1) = max(combined(:));
    else
        vec_of_maxima(1) = combined(lmindices(combined, 1, subject_mask));
    end
else
    vec_of_maxima = zeros(1, niters);
    start = 1;
end

random_berns = 2*(binornd(1,0.5, nsubj, niters )-1/2);
for iter = start:niters
    %     iter
    iter
    if mod(iter, 100) == 0
        iter
    end
    random_berns_for_iter = random_berns(:, iter);
    random_sample = find(random_berns_for_iter < 0);
    data2 = data;
    data2(:, random_sample) = -data(:, random_sample); %This doesn't subtract the mean though it probably should!
    
    if strcmp(stat, 'Z')
        combined = mean(data2, (D+1));
    elseif strcmp(stat, 'T')
        combined = mvtstat( data2 );
    end
    
    if isnan(subject_mask)
        vec_of_maxima(iter) = max(combined(:));
    else
        vec_of_maxima(iter) = combined(lmindices(combined, 1, subject_mask));
    end
end

threshold = prctile(vec_of_maxima, 100*(1-alpha) );

end