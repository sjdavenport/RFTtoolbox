function F = IntAlpha( voxmfd, type )
% Alpha( voxmfd, type ) computes a field .
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  voxmfd  an object of class VoxManifold
%  type    type of edge considered 'x', 'y' or 'z'
%--------------------------------------------------------------------------
% OUTPUT
%   F   
%--------------------------------------------------------------------------
% DEVELOPER TODOs:
%   - fix local isotropy assumption in L1 by computing the second integral
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow
%--------------------------------------------------------------------------

%% Check mandatory input and get important constants
%--------------------------------------------------------------------------

%% Main function
%--------------------------------------------------------------------------

switch type
    case 'x'
        
    case 'y'
        
    case 'z'
        
end

F = v.' * ( voxmfd.g * w );
           
return