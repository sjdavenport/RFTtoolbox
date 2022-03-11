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
mask_mult = 1 * mask;
if val ~= 0
    mask_mult( ~mask ) = val;
end

if val == -Inf || val == Inf
    % +0.2 is added to remove a bug, if val = -Inf/Inf and the outside of the
    % mask is zero!
    obj.field = repmat( mask_mult, [ ones( [ 1 D ] ), obj.fibersize ] ) .* (obj.field + 0.2) - 0.2;
else
    obj.field = repmat( mask_mult, [ ones( [ 1 D ] ), obj.fibersize ] ) .* obj.field;
end

% Put the mask into the object
obj.mask = mask;

return