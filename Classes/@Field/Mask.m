function obj = Mask( obj, val, mask )
% Mask( obj, mask ) constructs a Field object having white
% noise in the fiber.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  obj   of type Field
% Optional
%  val   a numeric used to fill the ~mask voxels
%  mask  a logical of size obj.sizeDomain.
%
%--------------------------------------------------------------------------
% OUTPUT
% obj  an object of class Fields representing white noise.
%
%--------------------------------------------------------------------------
% EXAMPLES
% sDomain = [ 10 10 ];
% mask = true( [ 10 10 ] );
% mask(:,2:4) = 0;
% lat_data = WhiteNoiseField( sDomain, 1, mask )
% 
% figure, clf,
% imagesc( lat_data.field ), colorbar
% 
% masked_lat_data = Mask( lat_data, mask )
% figure, clf,
% imagesc( masked_lat_data.field ), colorbar
% 
%--------------------------------------------------------------------------
% Author: Fabian Telschow
%--------------------------------------------------------------------------

%% Check optional input
%--------------------------------------------------------------------------

if ~exist( 'val', 'var' )
    val = 0;
end

if ~exist( 'mask', 'var' )
    mask = obj.mask;
end

% Check whether the mask and the field fit together
sMask  = size( mask );

if ~all( sMask == obj.masksize )
     error( "'mask' must be a logical array of size of the first length mask coordinates." )
end

%% Main function
%--------------------------------------------------------------------------
% Get dimension of the Domain
D = obj.D;

% Mask the field
mask_mult = repmat( mask, [ ones( [ 1 D ] ), obj.fibersize ] );

obj.field( ~mask_mult ) = val;

% Put the mask into the object
obj.mask = mask;

return