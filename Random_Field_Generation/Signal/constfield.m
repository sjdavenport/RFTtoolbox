function obj = constfield( A, mask_or_dim, xvals )
% constfield( A, mask_or_dim, xvals ) creates a Field object with constant
% value in the fiber.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  A             an F_1 x ... x F_K array
%  mask_or_dim   either a logical T_1 x ... x T_D array defining the mask
%                or a vector giving the size of the domain. In the latter
%                case the mask is given by true( mask_or_dim ).
% Optional
%  xvals         a 1 x D cell array containing the xvalues of the voxels.
%
%--------------------------------------------------------------------------
% OUTPUT
% obj   an object of class Field having constant values A in the fiber.
%
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% Author: Fabian Telschow
%--------------------------------------------------------------------------

%% Check mandatory input
%--------------------------------------------------------------------------

if ~isnumeric( A )
    error( "The input A needs to be a numerical array." )  
end

%% Check optional input
%--------------------------------------------------------------------------


%% Main function
%--------------------------------------------------------------------------

% Get size of input array
sA = size( A );
lA = length( sA );

% Create output field
obj = Field( mask_or_dim );

% Fill xvals, if provided
if exist( 'xvals', 'var' )
    obj.xvals = xvals;
end

% Fill the fiber by the array A
obj.field = permute( ...
                       repmat( A, [ ones( [ 1 lA ] ) obj.masksize ] ),...
                    [ ( lA + 1 ):(obj.D + lA), 1:lA ] );
return