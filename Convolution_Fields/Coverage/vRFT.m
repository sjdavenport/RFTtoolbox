 function [ output_image, threshold, max_finelat, L ] = vRFT(...
                            lat_data, params, alpha, version )
% vRFT_Field( lat_data, params, alpha, version ) runs voxelwise RFT 
% inference on a set of images to detect areas of activation using a 
% one-sample t-statistic.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  lat_data  an object of class Field containing the lattice data and the
%            mask
%  params    an object of class ConvFieldParams.
% Optional
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
%  output_image  the (fine lattice) output image
%  threshold     the voxelwise RFT threshold
%  max_finelat   the maximum on a fine lattice given by spacing
%--------------------------------------------------------------------------
% DEVELOPERS TODOS:
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport, Fabian Telschow
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Obtain the size of the lattice input
s_lat = size( lat_data );

% Obtain the dimensions of the data
Dim = s_lat( 1:end-1 );

% Obtain the number of dimensions
D = length( Dim );

%%  Add/check optional values
%--------------------------------------------------------------------------

if ~exist( 'alpha', 'var' )
    alpha = 0.05;
end

if ~exist( 'version', 'var' )
    switch lat_data.D
        case 1
            version = 0;
        case 2
            version = true;
        case 3
            version = [ true, true, true ];
    end
end

if ~exist( 'type', 'var' )
   % Default option of type
   type = "T";
end

if ~exist( 'df', 'var' )
   % Default option of type
   df = lat_data.fibersize-1;
end

%%  Main function
%--------------------------------------------------------------------------
% Calculate the convolution t field 
[ tfield_fine, cfields ] = convfield_t( lat_data, params );

% Estimate the LKCs of the convoltution field
dcfields = convfield( lat_data, params, 1 );
d2cfields = Field();

if lat_data.D == 3
    if version(3) == 1
        d2cfields = convfield( lat_data, params, 2 );
        [ L, L0 ] = LKC_voxmfd_est( cfields, dcfields, d2cfields,...
                                    version );
    else
        [ L, L0 ] = LKC_voxmfd_est( cfields, dcfields, d2cfields,...
                                    version );
    end
else
	[ L, L0 ] = LKC_voxmfd_est( cfields, dcfields, d2cfields, version );
end

% Calculate the threshold using the EEC heuristic
threshold = EECthreshold( alpha, L, L0, type, df );

% Determine the areas of the image where the t-field exceeds the threshold,
% need to use Newton Rhapson here
output_image = tfield_fine;

% Convert the output to double (can't remember why I do this!)
output_image.field = double( tfield_fine.field > threshold );

% add output of maximum, if you think we need it
max_finelat = 1;

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
