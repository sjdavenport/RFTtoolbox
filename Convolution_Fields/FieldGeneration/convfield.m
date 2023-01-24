function [ cfield, ss ] = convfield( lat_data, params, derivtype )
% CONVFIELD( lat_data, params, derivtype )
% generates an object of class Field containing the convolution field
% derived from lattice data smoothed with a seperable kernel.
% The generated field is evaluated on an equidistant grid with resolution
% increased by adding resadd voxels inbetween each voxel in each dimension.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  lat_data  an object of class Field with fiberD == 1.
%  params    an object of class ConvFieldParams
% Optional
%  derivtype  0/1/2, 0 returns the convolution field, 1 its derivative
%             and 2 its second derivative (at all points). Default is 0
%             i.e to return the field. Note that if D > 1, then derivtype = 2
%             has not yet been implemented.
%--------------------------------------------------------------------------
% OUTPUT
%  cfield   an object of class ConvField representing the field obtained by
%           smoothing the lat_data field object by the Kernel object.
%  ss       numeric to multiply the field by to make it variance 1.
%--------------------------------------------------------------------------
% EXAMPLES
% %% 1D 
% nvox = 10; FWHM = 2;
% lat_data = wfield(nvox); resadd = 0; enlarge = 0;
% params = ConvFieldParams( FWHM, resadd, enlarge);
% latcfield = convfield(lat_data, params);
% plot(latcfield.xvals{1}, latcfield.field, 'o-')
% hold on
% params.resadd = 10;
% cfield = convfield( lat_data, params);
% plot(cfield.xvals{1},cfield.field)
%
% %% 2D
% Dim = [10,10]; FWHM_x = 5; FWHM_y = 20;
% lat_data = wfield(Dim); resadd = 5; enlarge = 0;
% params = ConvFieldParams( [FWHM_x, FWHM_y], resadd, enlarge);
% latcfield = convfield(lat_data, params);
% imagesc(latcfield.field)
%
% Dim = [100,100]; FWHM = 5;
% lat_data = wfield(Dim); resadd = 1; enlarge = 0;
% params = ConvFieldParams( [FWHM, FWHM], resadd, enlarge);
% [latcfield, ss] = convfield(lat_data, params);
% latcfield.field = latcfield.field/sqrt(ss);
% imagesc(latcfield.field)
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow, Samuel Davenport                                            
%--------------------------------------------------------------------------

%% Check mandatory input
%--------------------------------------------------------------------------

% Allow for non field input
if ~isa( lat_data, 'Field' ) && isnumeric(lat_data)
    D = mask_dim(lat_data);
    if D == 1 && size(lat_data, 1) == 1
        lat_data = lat_data';
        lat_data = Field(lat_data, true(size(lat_data,1),1) );
    else
        lat_data = Field(lat_data, true(size(lat_data)));
    end
%     warning('PRobably should not enable this as can lead to errors?')
%     error('Need to sort the mask in one D here!')
end

% Allow for default FWHM input
if isnumeric(params)
    params = ConvFieldParams( repmat(params, 1, lat_data.D), 0 );
end

% Ensure that the kernel has the same dimension as the data itself, if not
% return an error
if length(params.kernel.kernel) ~= lat_data.D
    error('The dimensions of lat_data and params are not compatible')
end

% Check the lat_data input
if ~iscompletefield( lat_data )
    error( "lat_data must be a complete object of class Field." )
elseif lat_data.fiberD ~=1
    error( "The fiber dimension needs to be 1." )
else
    % Get dimension of the input data
    D = lat_data.D;
    % Get size of the mask of input data
    Dim = lat_data.masksize;
    % Get number of subjects
    nsubj  = lat_data.fibersize;
    % Get the xvals
    xvals  = lat_data.xvals;
    % Index for multidimensional coding
    indexD = repmat( {':'}, 1, D );
end

% Save the parameter from params in a small file
Kernel     = params.kernel;
resadd     = params.resadd;
enlarge    = params.enlarge;
lat_masked = params.lat_masked;

% Get the difference between voxels with resolution increase
dx = NaN * ones( [ 1 D ] );
for d = 1:D
    dx(d) = xvals{d}(2)-xvals{d}(1);
end

dx_hr = dx ./ ( resadd + 1 );
% These lines are pointless i.e. they do exactly the same thing as before
% small bug id resadd = 0 hence the if!
% if resadd ~= 0
%     dx_hr = dx ./ ( resadd + 1 );
% else
%     dx_hr = dx;
% end

% Reject input, if resadd is to large in 3D
if D == 3 && ( resadd > 18 )
    error( 'In 3D you shouldn''t use such high resolution for memory reasons' )
end

%% Add/check optional values
%--------------------------------------------------------------------------

%%% derivtype input
if ~exist( 'derivtype', 'var' )
    derivtype = 0;
end

%% Main function
%--------------------------------------------------------------------------

%%% Mask the lat_data if necessary
if lat_masked && ~lat_data.masked
    lat_data = Mask( lat_data );
end

%%% Get resolution increased size of the domain
% Dimensions for domain of the field with increased resolution
Dimhr = ( Dim - 1 ) * resadd + Dim; % Dimhr = Dim with high resolution

% Modify Dimhr by adding enlarge voxels to all sides
if enlarge ~= 0 
    if D == 1
        Dimhr = Dimhr + 2 * enlarge * [ 1, 0 ];
    else
        Dimhr = Dimhr + 2 * enlarge;
    end
end

%%% Expand lattice to new high resolution by filling in zeros
% Create index to fill the original data at the correct voxels of the high
% resolution data array
index = cell( [ 1 D ] );
for d = 1:D
    index{d} = ( enlarge + 1 ):( resadd + 1 ):( Dimhr(d) - enlarge );
end

% Increase the resolution of the raw data by introducing zeros
expanded_lat_data = squeeze( zeros( [ Dimhr, nsubj ] ) );
expanded_lat_data( index{:}, : ) = lat_data.field;

%%% Prepare output of class Field
% Create a ConvField class object with high resolution mask
cfield = ConvField();
% Fill the properties
cfield.kernel    = Kernel;
cfield.enlarge   = enlarge;
cfield.resadd    = resadd;
cfield.derivtype = derivtype;

cfield.mask      = mask_highres( lat_data.mask, resadd, enlarge );
xvals2 = cell( [ 1, D ] );

for d = 1:D
    xvals2{d} = ( xvals{d}(1) - enlarge * dx_hr(d) ) : dx_hr(d) : ...
                                    ( xvals{d}(end) + enlarge * dx_hr(d) );
end

cfield.xvals = xvals2; % This works just fine - the other option results in a bug
% if resadd ~= 0
%     cfield.xvals = xvals2;
% else
%     cfield.xvals = xvals;
% end

%%% Main loop: calculation of convolution fields
if D < 4
    % run the smoothing using fconv
    if derivtype == 0
        % calculates the convolution field
        [ tmp, ss ] = fconv( expanded_lat_data, Kernel.kernel, D,...
                              Kernel.truncation, dx_hr, Kernel.adjust );
        cfield.field = tmp;
        
    elseif derivtype == 1        
        % preallocate the output field for speed
        smooth_data = ones( [ Dimhr(1:D), nsubj, D ] );
        
        % Get the gradient object of the Kernel
        dKernel = Gradient( Kernel );
        
        % Calculates the derivatives of the convolution field
        for d = 1:D
            smooth_data(indexD{:},:,d) = fconv( expanded_lat_data,...
                                                 dKernel{d}.kernel, D, ...
                                                 dKernel{d}.truncation,...
                                                 dx_hr, dKernel{d}.adjust );
        end
        
        % Fill the derivatives into the cfield output
        cfield.field = smooth_data;
        
    elseif derivtype == 2
        % preallocate the output field for speed
        smooth_data = ones( [ Dimhr(1:D), nsubj, D, D ] );
        
        % Get the Hessian object of the Kernel
        Kernel = Hessian( Kernel );
        
        % Calculates all second derivatives of the convolution field
        for d = 1:D
           for dd = 1:D
                smooth_data(indexD{:},:,d,dd) = fconv( expanded_lat_data,...
                                                 Kernel{d,dd}.kernel, D, ...
                                                 Kernel{d,dd}.truncation,...
                                                 dx_hr, Kernel{d,dd}.adjust );
           end
        end
        
        % Fill the derivatives into the cfield output
        cfield.field = smooth_data;
        
    else
        error( 'Higher derivatives than 2 are not supported' )  
    end
    
else
    error('D != 1,2,3 has not been implemented yet!')
end

% @Fabian, I don't think this code makes sense, ss should be the same here
% as in fconv as ss stands for sum of squares. 
% if exist( 'ss', 'var' )
%     ss = 1 / sqrt(ss);
% else
%     ss = 1;
% end
% changes to:
if ~exist( 'ss', 'var' )
    ss = 1;
end

return