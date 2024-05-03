function [ output_image, threshold, maximum, L ] = vRFT_gauss(...
                                                  lat_data, params, alpha )
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
%  do2sample    0/1 whether or not to record the minimum as well as the
%               maximum in order to do a two sample test. Default is 0,
%               i.e. to do a one-sample test.
%  version   a logical/ logical vector. Length depends on voxmfd.D
%            - D = 1, always true.
%            - D = 2, either 1 or 0. 1 if L1 should be estimated, 0 else.
%            - D = 3, logical of length 3. version(1), indicates whether L2
%            should be estimated, version(2) whether the first integral is
%            used in L1 and version(3) whether the second integral is used.
%  alpha     the alpha level at which to threshold. Default is 0.05.
%            Recommend alpha <= 0.05 for best performance.
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


if ~exist( 'type', 'var' )
    % Default option of type
    field_type = "T";
end

if ~exist( 'df', 'var' )
    % Default option of type
    df = lat_data.fibersize-1;
end

% Set the default value of do2sample
if ~exist( 'do2sample', 'var')
    do2sample = 0;
end


%%  Main function
%--------------------------------------------------------------------------

%% LKC calculations
% Calculate the convolution t field
[ ~, cfields ] = convfield_t( lat_data, params );
G_cfields = Gaussianize(cfields);
tfield_fine = convfield_t(G_cfields, 0);

FWHM_est = mean(est_smooth(G_cfields.field)/(params.resadd + 1));
Lambda = FWHM2Lambda(FWHM_est);
Lambda_field = constfield( Lambda, G_cfields.masksize );
Lambda_field.mask = G_cfields.mask;
voxmfd  = VoxManifold( Lambda_field );
[ LKCs, L0 ] = LKC_est( voxmfd );

% Calculate the threshold using the EEC heuristic
L.L = real(LKCs); %Necessary as sometime there are small imaginary parts due to numerical inaccuracies
L.L0 = L0;
threshold = EECthreshold( alpha, L.L, L0, field_type, df );

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
