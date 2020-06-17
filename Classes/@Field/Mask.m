function obj = Mask( obj, mask )
% WhiteNoiseField( obj, FWHM ) constructs a Fields object having white
% noise in the fiber.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  obj   of type Field
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

%% Check mandatory input
%--------------------------------------------------------------------------

if ~all( islogical( mask(:) ) )
     error( "'mask' must be a logical array." )
end

% Check whether the mask and the field fit together
sMask  = size( mask );
sField = size( obj.field );
if ~all( sMask == sField( 1:length( sMask ) ) )
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
obj.field = repmat( mask, [ ones( [ 1 D ] ), sField( D+1:end ) ] ) .* obj.field;

% Put the mask into the object
obj.mask = mask;

return