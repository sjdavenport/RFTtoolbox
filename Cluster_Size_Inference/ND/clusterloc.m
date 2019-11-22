function clusters = clusterloc(data, thresh, connectivity_criterion)
% clusterloc(data, thresh, connectivity_criterion) calculates the number
% of connected components in an array that lie above a threshold.
%--------------------------------------------------------------------------
% REQUIRES the image processing toolbox for Matlab, specifically the
% function bwconncomp.m .
%--------------------------------------------------------------------------
% ARGUMENTS
% data      a 2 or 3 dimensional array of real values.
% thresh    a threshold
% connectivity_criterion    To Add!
%--------------------------------------------------------------------------
% OUTPUT
% number_of_clusters    The number of connected components above the
%                       threshold thresh
% occurences
% sizes
%--------------------------------------------------------------------------
% EXAMPLES
% %2D
% thresh = 1;
% sims = randn(10,10)
% clusters = clusterloc(sims, thresh);
% sims > thresh
% clusters
% 
% %3D
% thresh = 1;
% sims = randn(10,10,10)
% clusters = clusterloc(sims, thresh);
% sims > thresh
% clusters
%--------------------------------------------------------------------------
% AUTHORS: Samuel Davenport
s = size(data);
ones_and_zeros = zeros(s);

if s(1) == 1
    D = 1;
else
    D = length(s);
end
ones_and_zeros(data > thresh) = 1;

%If we're in the 2d case only look at components connected vertically
%and horizontally (i.e. no diagonals!)
if D == 2
    if nargin < 3
        connectivity_criterion = 4;
    end
    if sum(connectivity_criterion == [4,8]) ~= 1
        error('In 2D the connectivity criterion must be 4 or 8')
    end
    conComponents = bwconncomp(ones_and_zeros, connectivity_criterion);
elseif D == 3
    if nargin < 3
        connectivity_criterion = 18;
    end
    if sum(connectivity_criterion == [6,18,26]) ~= 1
        error('In 3D the connectivity criterion must be 6, 18 or 26')
    end
    conComponents = bwconncomp(ones_and_zeros, connectivity_criterion);
else
    error('The dimension must be 2 or 3.');
end

%Return the true number of connected components/clusters.
clusters.nC = conComponents.NumObjects;
clusters.index_locations = conComponents.PixelIdxList;
clusters.cextent = cellfun(@(x) numel(x), clusters.index_locations);

clusters.edgeind = zeros(1, clusters.nC);
clusters.arraylocs = cell(1, clusters.nC);
for c = 1:clusters.nC
    clocs = clusters.index_locations{c};
    clusters.arraylocs{c} = zeros(length(clocs), D);
    if D == 2
        for I = 1:length(clocs)
            [clusters.arraylocs{c}(I, 1),clusters.arraylocs{c}(I, 2)] = ind2sub(s, clocs(I));
        end
        if any(clusters.arraylocs{c}(:, 1) == 1) || any(clusters.arraylocs{c}(:, 2) == 1) || any(clusters.arraylocs{c}(:, 1) == s(1)) || any(clusters.arraylocs{c}(:, 2) == s(2))
            clusters.edgeind(c) = 1;
        else
            clusters.edgeind(c) = 0;
        end
    elseif D == 3
        for I = 1:length(clocs)
            [clusters.arraylocs{c}(I, 1),clusters.arraylocs{c}(I, 2), clusters.arraylocs{c}(I, 3)] = ind2sub(s, clocs(I));
        end
        if any(clusters.arraylocs{c}(:, 1) == 1) || any(clusters.arraylocs{c}(:, 2) == 1) || any(clusters.arraylocs{c}(:, 3) == 1) || any(clusters.arraylocs{c}(:, 1) == s(1)) || any(clusters.arraylocs{c}(:, 2) == s(2)) || any(clusters.arraylocs{c}(:, 3) == s(3))
            clusters.edgeind(c) = 1;
        else
            clusters.edgeind(c) = 0;
        end
    end
end

end

% clusters.edgeind(c) = 0;
% for I = 1:length(clocs)
%     [x,y,z] = ind2sub(s, clocs(I));
%     if x == 1 || y == 1 || z == 1 || x == s(1) || y == s(2) || z == s(3)
%         clusters.edgeind(c) = 1;
%         break
%     end
% end