function [lat_data, standardized_field, standard_data] = ...
                                  Gaussianize( lat_data, stdsmo, usetrans )
% GAUSSIANIZE( lat_data ) takes in a field of data on a lattice and 
% Gaussianizes it. To do this it takes the data and demeans and 
% standardizes it voxelwise and combines over voxels to obtain a null
% marginal distribution. It then compares the original, standardized 
% but not demeaned, data to this null marginal distribution calculates a 
% quantile and Gaussianizes it via comparison to the Gaussian pdf. 
% The aim of this is to aim Gaussian distributions at each voxel (note 
% however that the field transformed this way are only approximately Gaussian
% given a finite number of subjects.). One of the benefits of this 
% transformation is to ensure that the tails do not overly influence the data.
% This is especially important in fMRI where due to artefacts in the data 
% it is common that some voxels of some subjects have very high tails so 
% that an indivdiual subject can dominate the distribution
% at a given voxel. This transformation increases the power to detect
% activation and ensures that the assumptions of multiple testing methods
% hold. Note that it is helpful for both the assumptions of permutation and 
% random field theory.
%
% In small samples the data is only approximately Gaussian as
% (X-muhat)/sigmahat has a slightly different distribution from X/sigmahat.
% This is typically only an issue in data that has an extreme amount of
% weight near zero as then values of muhat are highly variable. 
%
% If test-statistics are being used then Gaussianity can further be ensured
% by the central limit theorem. See Davenport (2021) for details.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  lat_data    a field of data on a lattice. For best performance this
%              field should be made up of at least 10 subjects. 
% Optional
%  stdsmo      a smoothing parameter for the standard deviation. Default is
%              not to smooth, i.e. to set stdsmo = 0.
%  usetrans    different transformation options (to be specified!)
%--------------------------------------------------------------------------
% OUTPUT
%  lat_data    the Gaussianized field of data
%  standardized_field  (X_n - muhat)/sigmahat for 1 <= n <=leq N
%  standard_data  X_n/sigmahat for 1 <= n <= N
%--------------------------------------------------------------------------
% EXAMPLES
% lat_data = wfield( [20,20], 100, 'T', 3 )
% gaussianized_data = Gaussianize( lat_data )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist('stdsmo', 'var')
    stdsmo = 0;
end

if ~exist('usetrans', 'var')
    usetrans = 0;
end

%% Allow for non field input
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

% Ensure that the field is masked
lat_data = Mask(lat_data);

%%  Main Function Loop
%--------------------------------------------------------------------------

% Standardize
% std_dev = std(lat_data.field, 0, lat_data.D + 1);
std_dev = std(lat_data);

% Smooth the standard deviation if that is included as an option
% Need to work on this and explore it further it doesn't work atm
if stdsmo > 0 
    params = ConvFieldParams(repmat(stdsmo, 1 , lat_data.D), 0, 0);
    std_dev = Mask(convfield(std_dev, params));
    onefield = std_dev./std_dev.field;
    twofield = Mask(convfield(onefield, params));
    std_dev = std_dev./twofield;
    error('Weird things happen when you try to smooth the MNI mask, sum(nonanlocs(:)) is way too small!')
end

% Calculate the mean at each voxel
mean_dev = mean(lat_data.field, lat_data.D + 1);

if stdsmo >= 0
    % Standardize without demeaning
    standard_data = lat_data.field./std_dev.field;
    % Standardize after demeaing to get the null distribution
    standardized_field = (lat_data.field-mean_dev)./std_dev.field;
else
    % Obtain orig data
    standard_data = lat_data.field;
    % Demean
    standardized_field = lat_data.field-mean_dev;
end

% Obtain the locations of the data that lies within the mask
nonnanlocs = ~isnan(standard_data);

if usetrans == 1
    % Work in progress: this loop corrects for the fact that you have 
    % subtracted the mean and ensures that normal input gives normal output.
    if lat_data.fibersize == 1
        return
    else
        data_vc = vec_ecdf( standard_data(nonnanlocs), standardized_field(nonnanlocs) );
        load(['C:\Users\12SDa\davenpor\Data\RestingStateData\EklundData\Gaussianize_transmutes\vc2_nsubj_', ...
            num2str(lat_data.fibersize)], 'vc_dist')
    end
    normalized_data = norminv(vec_ecdf( data_vc, vc_dist ), 0, 1);
elseif usetrans == 0
    % Default loop
    % vec_ecdf calculates the empirical cdf of the standard_data based
    % on the standardized_data then norminv is applied to convert the
    % quantiles to normal ones and thus Gaussianize the data
    normalized_data = norminv(vec_ecdf( standard_data(nonnanlocs), standardized_field(nonnanlocs) ), 0, 1);
elseif usetrans == 2
    % Use the square root transform rather than the estimated one. The
    % utility of this and other transformations needs to be investigated as
    % it gives good performance despite being quite simple.
    normalized_data = abs(standard_data(nonnanlocs)).^(1/2).*sign(standard_data(nonnanlocs));
end

% Set the data within the mask to be the Gaussianized data
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
