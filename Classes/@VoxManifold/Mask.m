function obj = Mask( obj, mask )
% Mask( obj, mask ) masks all the Field objects in the VoxelManifold object.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  obj   of type Field
% Optional
%  mask  a logical of size obj.sizeDomain.
%
%--------------------------------------------------------------------------
% OUTPUT
% obj  an object of class Fields representing white noise.
%
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% Author: Fabian Telschow
%--------------------------------------------------------------------------

%% Check optional input
%--------------------------------------------------------------------------

if ~exist( 'mask', 'var' )
    mask = obj.mask;
end

% Check whether the mask and the field fit together
sMask  = size( mask );
sField = size( obj.mask );
if ~all( sMask == sField )
     error( "'mask' must be a logical array of size of the first length mask coordinates." )
end

%% Main function
%--------------------------------------------------------------------------
% Get dimension of the Domain
D = length( sMask );
if any( sMask == 1) && D == 2
    D = 1;
end

% Mask the field
obj.g.field = repmat( mask, [ ones( [ 1 D ] ), sField( D+1:end ) ] )...
                                                            .* obj.g.field;

% Put the mask into the object
obj.mask = mask;
obj.g.mask = mask;

return