function dx = get_dx( obj )
% get_dx( obj ) obtains a vector giving the dx in each direction of the grid.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  obj   of type Field
%
%--------------------------------------------------------------------------
% OUTPUT
% dx  1 x obj.D vector containing the dx in each dimension.
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
dx = ones( [ 1 obj.D ] );
for d = 1:obj.D
    dx(d) = obj.xvals{d}(2)-obj.xvals{d}(1);
end

return