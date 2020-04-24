function [ out ] = perm_cluster( data, stat, niters, include_original, subject_mask, alpha )
% perm_cluster uses the permutation test to perform clustersize inference.
%--------------------------------------------------------------------------
% ARGUMENTS
% 
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Sam Davenport.
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

[nsubj, ~] = size(data);

if include_original
%     niters = niters - 1;
    start = 2;
    
    if strcmp(stat, 'Z')
        combined = mean(data);
    elseif strcmp(stat, 'T')
        combined = mvtstat( data);
    end
    
    max_cluster_sizes = zeros(1, niters );
    if isnan(subject_mask)
        max_cluster_sizes(1) = max(combined(:));
    else
        max_cluster_sizes(1) = combined(lmindices(combined, 1, subject_mask));
    end
else
    max_cluster_sizes = zeros(1, niters);
    start = 1;
end

random_berns = 2*(binornd(1,0.5, nsubj, niters )-1/2);

for iter = start:niters
%     iter
    iter
    random_berns_for_iter = random_berns(:, iter);
    random_sample = find(random_berns_for_iter < 0);
    data2 = data;
    data2(random_sample, :) = -data(random_sample, :); %This doesn't subtract the mean though it probably should!
    
    if strcmp(stat, 'Z')
        combined = mean(data2);
    elseif strcmp(stat, 'T')
        combined = mvtstat( data2 );
    end
    
    if isnan(subject_mask)
        max_cluster_sizes(iter) = max(combined(:));
    else
        max_cluster_sizes(iter) = combined(lmindices(combined, 1, subject_mask));
    end
end

threshold = prctile(max_cluster_sizes, 100*(1-alpha) );

end

