function lat_data = Gaussianize( lat_data )
% GAUSSIANIZE( lat_data ) takes in a field of data on a lattice and 
% Gaussianizes it demeans. To do this it takes the data and demeans and 
% standardizes it voxelwise and combines over voxels to obtain a null
% marginal distribution (it makes the assumption that the marginal 
% distribution a condition that is of course much weaker than Gaussianity).
% It then compares the original, standardized but not demeaned, data
% to this null marginal distribution calculates a quantile and Gaussianizes
% it via comparison to the Gaussian pdf. The effect of this is to ensure
% marginal Gaussian distributions at each voxel (note however that the
% field transformed this way are not necessarily Gaussian). One of the 
% benefits of this transformation is to ensure that the tails do not overly
% influence the data. This is especially important in fMRI where due to
% artefacts in the data it is common that some voxels of some subjects have 
% very high tails so that an indivdiual subject can dominate the distribution
% at a given voxel. This transformation increases the power to detect
% activation and ensures that the assumptions of multiple testing methods
% hold. Note that it is helpful for both the assumptions of permutation and 
% random field theory.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  lat_data    a field of data on a lattice. For best performance this
%              field should be made up of at least 10 subjects. 
%--------------------------------------------------------------------------
% OUTPUT
%  lat_data    the Gaussianized field of data
%--------------------------------------------------------------------------
% EXAMPLES
% lat_data = wfield( [20,20], 100, 'T', 3 )
% gaussianized_data = Gaussianize( lat_data )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

% Allow for non field input
if ~isa( lat_data, 'Field' ) && isnumeric(lat_data)
    temp_lat_data = lat_data;
    s_lat_data = size(lat_data);
    s_lat_data = s_lat_data(1:end-1);
    if length(s_lat_data) == 1
        s_lat_data = [s_lat_data, 1];
    end
    lat_data = Field(true(s_lat_data));
    lat_data.field = temp_lat_data;
    clear temp_lat_data;
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
