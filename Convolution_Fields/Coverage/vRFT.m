function [ output_image, threshold, maximum, L ] = vRFT(...
                        lat_data, params, ninitpeaks, version, alpha, L0 )
% VRFT( lat_data, params, L0, ninitpeaks, alpha, version ) runs voxelwise
% RFT inference on a set of images to detect areas of activation using a
% one-sample t-statistic.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  lat_data  an object of class Field containing the lattice data and the
%            mask
%  params    an object of class ConvFieldParams.
% Optional
%  ninitpeaks   the number of initial peaks to use to find the maximum of
%               the convolution field
%  alpha     the alpha level at which to threshold. Default is 0.05.
%            Recommend alpha <= 0.05 for best performance.
%  version   a logical/ logical vector. Length depends on voxmfd.D
%            - D = 1, always true.
%            - D = 2, either 1 or 0. 1 if L1 should be estimated, 0 else.
%            - D = 3, logical of length 3. version(1), indicates whether L2
%            should be estimated, version(2) whether the first integral is
%            used in L1 and version(3) whether the second integral is used.
%--------------------------------------------------------------------------
% OUTPUT
%  output_image   a logical array (with dimensions that of the resolution
%                 increased data) with 0 if a given point is below the
%                 threshold and 1 if above
%  threshold     the voxelwise RFT threshold
%  maximum    a structure with fields:
%   maximum.lat      the maximum on the original lattice
%   maximum.finelat  the maximum on the resolution increased lattice
%   maximum.conv     the maximum of the convolution field calculated using
%                   an optimization algorithm (findconvpeaks)
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport, Fabian Telschow
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Get the number of dimensions
D = lat_data.D;

%%  Add/check optional values
%--------------------------------------------------------------------------

if ~exist( 'alpha', 'var' )
    alpha = 0.05;
end

% If this is not supplied it searches for the max using 1 initial peak
if ~exist( 'ninitpeaks', 'var' )
    ninitpeaks = 1;
end

if ~exist( 'version', 'var' )
    switch lat_data.D
        case 1
            version = false;
        case 2
            version = true;
        case 3
            version = [ true, true, false ];
            %             version = [ true, true, true ];
    end
end

if ~exist( 'type', 'var' )
    % Default option of type
    field_type = "T";
end

if ~exist( 'df', 'var' )
    % Default option of type
    df = lat_data.fibersize-1;
end

if ~exist( 'L0', 'var')
    L0 =  EulerChar(lat_data.mask, 0.5, lat_data.D);
end

%%  Main function
%--------------------------------------------------------------------------

%% LKC calculations
% Calculate the convolution t field
[ tfield_fine, cfields ] = convfield_t( lat_data, params );

if isnumeric(version)
    L = version;
elseif isstr(version) && strcmp(version, 'HPE')
    HPE  = LKC_HP_est( cfields, 1, 1 );
    L = HPE.hatL;
else
    % Estimate the LKCs of the convoltution field
    dcfields = convfield( lat_data, params, 1 );
    d2cfields = Field();
    
    % @Fabian, L0 is buggy as it's calculated on the resadd increased mask.
    % Plus this is extra computation time so perhaps by default it shouldn't be computed?
    % As it will be the same every time so it doesn't make sense to recompute
    % it
    if lat_data.D == 3
        if version(3) == 1
            d2cfields = convfield( lat_data, params, 2 );
            [ L, L0bug ] = LKC_voxmfd_est( cfields, dcfields, d2cfields,...
                version );
        else
            [ L, L0bug ] = LKC_voxmfd_est( cfields, dcfields, d2cfields,...
                version );
        end
    else
        [ L, L0bug ] = LKC_voxmfd_est( cfields, dcfields, d2cfields, version );
    end
end

% Calculate the threshold using the EEC heuristic
L = real(L);
threshold = EECthreshold( alpha, L, L0, field_type, df );

% Determine the areas of the image where the t-field exceeds the threshold,
% need to use Newton Rhapson here
output_image = tfield_fine;

% Convert the output to double (can't remember why I do this!)
output_image.field = double( tfield_fine.field > threshold );

%% Obtain the maximum
% Define the locations of the original voxels (in terms of the coordinates
% of the resolution increased values.
orig_lattice_locs = cell(1,D);
for d = 1:D
    orig_lattice_locs{d} = (ceil(params.resadd/2) + 1):(params.resadd+1):length(cfields.xvals{d});
end

% Calculate the maximum on the original lattice
tfield_lat = tfield_fine.field(orig_lattice_locs{:});
maximum.lat = max(tfield_lat(:).*zero2nan(lat_data.mask(:)));

% Calculate the maximum on the fine lattice
maximum.finelat = max(tfield_fine.field(:).*zero2nan(tfield_fine.mask(:)));

if ninitpeaks > 0
    % Find the locations of the top peaks of the fine tfield on the mask
    % (these are needed to initialize the peak finding algorithm)
    peak_est_locs = lmindices(tfield_fine.field, ninitpeaks, tfield_fine.mask);
    
    % Ensure that peak_est_locs is D by npeaks
    if D == 1
        peak_est_locs = peak_est_locs';
    end
    
    % Covert to the correct coordinate system (determined by tfield_fine.xvals)
    peak_est_locs = xvaleval(peak_est_locs, tfield_fine.xvals);
    
    % Convert to a cell array in 1D for input to findconvpeaks
    if D == 1
        peak_est_locs = num2cell(peak_est_locs);
    end
    
    % Obtain the FWHM from the kernel (need to update findconvpeaks to work
    % with an arbitrary kernel!)
    FWHM = Gker2FWHM( params.kernel );
    
    % Calculate the maximum of the convolution field
    [~, max_tfield_at_lms] = findconvpeaks(lat_data.field, FWHM, peak_est_locs, 'T', lat_data.mask);
    maximum.conv = max(max_tfield_at_lms);
    
    maximum.allmaxima = max_tfield_at_lms;
    
    % There are some residual mismatches in 1D between convfield and
    % applyconvfield that need to be resolved.
    if D == 1
        maximum.conv = max(maximum.conv, maximum.finelat);
    end
else
    maximum.conv = maximum.finelat;
    maximum.allmaxima = maximum.finelat;
end

end

% DEPRECATED
% high_local_maxima = lmindices(tfield_fine, 3);
% Calculate initial estimates of peak location
% if D == 1
%     peak_est_locs = [NaN,setdiff(xvals_fine(high_local_maxima),[1,nvox])];
% end

% tcf = @(x) tcfield( x, lat_data, FWHM );

% if length(peak_est_locs) == 1
%     tfield_at_lms = -Inf; %If the local max occurs at the boundary you don't need to account for it.
% else
%     top_lmlocs = findconvpeaks_t(lat_data, FWHM, peak_est_locs);
%     tfield_at_lms = tcf(top_lmlocs);
% end
%
% % Calculate the maximum on the lattice and of the convolution field
% max_finelat = max(tfield_fine);
% max_conv = max([tfield_at_lms,max_finelat]); %Included for stability in case the maximum finding didn't work correctly.
