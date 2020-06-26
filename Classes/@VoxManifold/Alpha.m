function F = Alpha( voxmfd, type )
% Alpha( voxmfd, type ) computes a field .
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  voxmfd  an object of class VoxManifold
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

F = v.' * ( voxmfd.g * w );
           
return