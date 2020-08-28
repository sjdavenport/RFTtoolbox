function lat_data = Gaussianize( lat_data )
% Gaussianize takes in data on a lattice 
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% lat_data = wtfield( [20,20], 100, 3 )
% gaussianized_data = Gaussianize( lat_data )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'opt1', 'var' )
   % default option of opt1
   opt1 = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------

% Standardize
std_dev = std(lat_data.field, 0, lat_data.D + 1);
mean_dev = mean(lat_data.field, lat_data.D + 1);

standard_data = lat_data.field./std_dev; 
nonnanlocs = ~isnan(standard_data);
standardized_field = (lat_data.field-mean_dev)./std_dev; 

normalized_data = norminv(vec_ecdf( standard_data(nonnanlocs), standardized_field(nonnanlocs) ), 0, 1);
lat_data.field(nonnanlocs) = normalized_data;

end

%Old -1
% [standard_data, sort_index] = sort(standard_data);
% sort_index = invPerm(sort_index);
% standard_data = standard_data(~isnan(standard_data));
% standard_data = [norminv(vec_ecdf( standard_data, dist_vec ), 0, 1), nan(1, prod(fieldsize) - n_nonnans)] ;
% 
% Gaussianized_data = Field(mask);
% Gaussianized_data.field = reshape(standard_data(sort_index),fieldsize);


% Old 0
% G_data = nan(1, length(sort_index));
% G_data(1:n_nonnans) = (0:(n_nonnans-1))/n_nonnans;
% G_data(1) = 1/n_nonnans;
% G_data = G_data(invPerm(sort_index));
% G_data = reshape(G_data, fieldsize);

% Old 1
% G_data = sum(vectorized_data < vectorized_data');

% Old 2
% G_data = zeros(1, length(vectorized_data));
% 
% for I = 1:length(vectorized_data)
%     G_data(I) = sum(vectorized_data <= vectorized_data(I));
% % end
% G_data = G_data/length(vectorized_data);
% G_data(G_data == 1) = 1-1/length(vectorized_data);

% Old 3
% indicator_sum = zeros(size(lat_data.field));
% 
% for I = 1:length(vectorized_data)
% %     modul(I,1000)
%     indicator_sum = indicator_sum + (lat_data.field  <= vectorized_data(I));
% end
% indicator_sum = indicator_sum/length(vectorized_data);
% indicator_sum(indicator_sum == 0) = 1/length(vectorized_data)/2;
% indicator_sum(indicator_sum == 1) = 1-1/length(vectorized_data)/2;
