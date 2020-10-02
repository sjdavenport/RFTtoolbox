function obj = cut2mask( obj )
% cut2mask( obj ) returns a Field class object, where the domain is cutted
% to the smalest box containing the mask.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  obj   object of class Field
%
%--------------------------------------------------------------------------
% OUTPUT
% cobj  object of class Field where the domain is cutted to the 
%
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% Author: Fabian Telschow
%--------------------------------------------------------------------------

%% Check optional input
%--------------------------------------------------------------------------

%% Main function
%--------------------------------------------------------------------------
% Get the mask
mask  = obj.mask;

% Find the defining box
vals = cell( [ 1 obj.D ] );

switch obj.D
    case 1
        tmp     = find( sum( mask ) ~= 0 );
        vals{1} = tmp( 1 ) : tmp( end );
        obj    = obj( vals{1} );
    case 2
        tmp     = find( sum(sum( mask, 3 ), 2) ~= 0 );
        vals{1} = tmp( 1 ):tmp( end );
        tmp     = find( sum(sum( mask, 3 ), 1) ~= 0 );
        vals{2} = tmp( 1 ):tmp( end );
        obj    = obj( vals{1}, vals{2} );
    case 3
        tmp     = find( sum(sum( mask, 3 ), 2) ~= 0 );
        vals{1} = tmp( 1 ):tmp( end );
        tmp     = find( sum(sum( mask, 3 ), 1) ~= 0 );
        vals{2} = tmp( 1 ):tmp( end );
        tmp     = find( sum(sum( mask, 2 ), 1) ~= 0 );
        vals{3} = tmp( 1 ):tmp( end );
        obj    = obj( vals{1}, vals{2}, vals{3} );
end
return