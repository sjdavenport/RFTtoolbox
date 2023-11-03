function [ surviving_cluster_im, surviving_clusters, surviving_clusters_vec] = ...
                             cluster_im( dim, index_locations, threshold )
% CLUSTER_IM( dim, index_locations, threshold )
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'opt1', 'var' )
   % Default value
   opt1 = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
surviving_cluster_im = zeros(dim);
surviving_clusters_vec = {};
nsurvivors = 0;
for I = 1:length(index_locations)
  if length(index_locations{I}) > threshold
      nsurvivors = nsurvivors + 1;
      surviving_cluster_im(index_locations{I}) = 1;
      surviving_clusters_vec{nsurvivors} = index_locations{I};
  end
end
surviving_clusters = convindall(surviving_clusters_vec);
end

