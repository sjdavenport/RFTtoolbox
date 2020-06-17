function obj = WhiteNoiseField( sizeDomain, sizeFiber, mask )
% WhiteNoiseField( obj, FWHM ) constructs a Fields object having white
% noise in the fiber.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  sizeDomain  a vector containing the size of the domain
%
% Optional
%  sizeFiber   a vector containing the size of the fiber. Default is 1,
%              i.e. the field is scalar.
%  mask        a logical of size sizeDomain. Default true( sizeDomain ).
%
%--------------------------------------------------------------------------
% OUTPUT
% obj  an object of class Fields representing white noise.
%
%--------------------------------------------------------------------------
% EXAMPLES
% %% % Simple example with whole domain mask
% %% scalar field
% lat_data = WhiteNoiseField( [4 2 3] )
%
% %% many subjects field
% lat_data = WhiteNoiseField( [4 2 3], 100 )
%
% %% Simple example with mask
% sDomain = [ 4 2 3 ]
% mask = true( [ 4 2 3 ] )
% mask(:,1,2) = 0;
% lat_data = WhiteNoiseField( [4 2 3], 1, mask )
% 
%--------------------------------------------------------------------------
% Author: Fabian Telschow
%--------------------------------------------------------------------------

%% Check mandatory input
%--------------------------------------------------------------------------
% Check input argument and make it a vector
if  ~isnumeric( sizeDomain )
    error( "The sizeDomain must be a 1xD or Dx1 numerical vector." )
else
    tmp = size( sizeDomain );
    if ~( length( tmp ) == 2 && any( tmp == 1 ) )
        error( "The sizeDomain must be a 1xD or Dx1 numerical vector." )
    elseif tmp( 2 ) == 1
        sizeDomain = sizeDomain';
    end
end

%% Check mandatory input
%--------------------------------------------------------------------------

% Fill the defaults of the optional parameters if necessary
if ~exist( 'sizeFiber', 'var')
    sizeFiber = 1;
end

if ~exist( 'mask', 'var')
    mask = true( sizeDomain );
end

% Check the optional inputs
if  ~isnumeric( sizeFiber )
    error( "The sizeFiber must be a 1xK or Kx1 numerical vector." )
else
    tmp = size( sizeFiber );
    if ~( length( tmp ) == 2 && any( tmp == 1 ) )
        error( "The sizeFiber must be a 1xK or Kx1 numerical vector." )
    elseif tmp( 2 ) == 1
        sizeFiber = sizeFiber';
    end
end

if ~all( islogical( mask(:) ) )
     error( "'mask' must be a logical array of size sizeDomain." )
end

if ~all( size( mask ) == sizeDomain )
     error( "'mask' must be a logical array of size sizeDomain." )
end

%% Main function
%--------------------------------------------------------------------------
obj = Field();
obj.field = randn( [ sizeDomain sizeFiber ] );
obj.mask = mask;

return