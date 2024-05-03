function lat_data = apower_trans( lat_data, power )
% ASINH_TRANS( lat_data ) transforms data with the inverse hyperbolic sinh
% transformation.
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
% gaussianized_data = asinh_trans( lat_data )
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

% Standardize without demeaning
standard_data = lat_data.field./std_dev.field; 

lat_data.field = apower(standard_data, power);
% lat_data.field = apower(standard_data);

end

