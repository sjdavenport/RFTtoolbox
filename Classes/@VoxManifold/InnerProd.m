function F = InnerProd( voxmfd, v, w )
% InnerProd( voxmfd, v, w ) computes the inner product of the vector fields
% v and w with respect to the Riemannian metric from voxmfd.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  voxmfd  an object of class VoxManifold
%--------------------------------------------------------------------------
% OUTPUT
%   F   an object of class field representing the scalar field obtained by
%       v.' * voxmfd.g * w
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