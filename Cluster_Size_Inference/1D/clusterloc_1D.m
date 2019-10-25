function clusters = clusterloc_1D(data, u, xvals)
% oneDclusterloc(data, u, xvals) find the clusters in a 1D random field.
%--------------------------------------------------------------------------
% ARGUMENTS
% data      a 1 dimensional array of real values
% u         a threshold
% xvals     a vector giving the location of each voxel in space
%--------------------------------------------------------------------------
% OUTPUT
% clusters a structure with fields:
%
% nC        the number of clusters above the threshold
% locs      an nC length cell array where each entry is a 2 point vector
%           with the start and end point of each cluster in terms of the
%           xvalues
% indices   an nC length cell array where each entry is a 2 point vector
%           giving the indicies of the start and end point of each cluster 
% maxloc    and nC length vector giving the location of the maximum within
%           each cluster.
% peakinclus    an nC length vector with entries 0 and 1, 1 meaning
%               that that the maximum of the cluster is a local maximum
%               (and is not at the boundary).
% extent    an nC length vector giving the cluster extent of each cluster
% edgeindicator     an nC length vector with entries 0 and 1, 1 meaning 
%                   that the cluster lies at the edge.
% data      the input data
% xvals     the input xvals    
%--------------------------------------------------------------------------
% EXAMPLES
% L = 10;
% data = normrnd(0,1,1,L);
% [cfield, xvals] = convfield_1D(data, 1);
% thresh = 1;
% plot(xvals, cfield)
% hold on
% plot(xvals, thresh*ones(1, length(xvals)), 'Linewidth', 2)
% clusters = clusterloc_1D(cfield, thresh)
%--------------------------------------------------------------------------
% AUTHOR: Samuel J. Davenport
if nargin < 3
    xvals = 1:length(data);
elseif length(xvals) ~= length(data)
    error('The length of xvals must be the same as the length of the data')
end
delta = xvals(2) - xvals(1);

above_u = data > u;

diff_above_u = diff(above_u);

start_locs = find(diff_above_u == 1)' + 1;
end_locs = find(diff_above_u == -1)';

edges = [0,0];
if above_u(1) == 1
    start_locs = [1,start_locs];
    edges(1) = 1;
end
if above_u(end) == 1
    end_locs = [end_locs, length(data)];
    edges(2) = 1;
end

clusters.nC = length(start_locs);
if ~isequal(clusters.nC, length(end_locs))
    error('for some reason the number of start locs is not the same as the number of end locs')
end
clusters.locs = cell(clusters.nC, 1);
clusters.indices = cell(clusters.nC, 1);
clusters.maxloc = zeros(1,clusters.nC);
clusters.peakinclus = zeros(1,clusters.nC);
clusters.cextent = zeros(1,clusters.nC);

for J = 1:clusters.nC
    clusters.indices{J} = [start_locs(J), end_locs(J)];
    clusters.locs{J} = xvals(clusters.indices{J}); %[xvals(start_locs(J)), xvals(end_locs(J))];
    cluster_indices = start_locs(J):end_locs(J);
    [~,maxidx] = max(data(cluster_indices)); 
    clusters.maxloc(J) = xvals(cluster_indices(maxidx));
    
    if (clusters.maxloc(J) > xvals(start_locs(J))) && (clusters.maxloc(J) < xvals(end_locs(J)))
        clusters.peakinclus(J) = 1;
    end
    clusters.cextent(J) = xvals(end_locs(J)) - xvals(start_locs(J)) + delta;
end

clusters.edgeindicator = zeros(1,clusters.nC);
if edges(1) == 1
    clusters.edgeindicator(1) = 1;
end
if edges(2) == 1
    clusters.edgeindicator(end) = 1;
end

clusters.data = data;
clusters.xvals = xvals;
end