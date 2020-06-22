function obj = EuclideanMetric( obj, mask )
% EuclideanMetric( obj ) fills the field property of a Field object such
% that it has the euclidean metric as fibers.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  obj   of type Field
% Optional
%  mask  a logical of size obj.sizeDomain.
%
%--------------------------------------------------------------------------
% OUTPUT
% obj  an object of class Fields representing the Euclidean metric.
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

%% Main function
%--------------------------------------------------------------------------
% Get dimension of the Domain
sMask = size( mask );
D = length( sMask );
if any( sMask == 1) && D == 2
    D = 1;
end

% Mask the field
obj.field = permute( repmat( eye( obj.D ), [ 1 1 sMask ] ), ...
                                                        [ 3:(D+2), 1, 2 ] );

% Put the mask into the object
obj.mask = mask;

return