function F = InnerProd( voxmfd, v, w )
% InnerProd( voxmfd, v, w ) computes the inner product of the vector fields
% v and w with respect to the Riemannian metric from voxmfd.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  voxmfd  an object of class VoxManifold
%  v       an object of class Field with fibersize = voxmfd.D
%  w       an object of class Field with fibersize = voxmfd.D
%--------------------------------------------------------------------------
% OUTPUT
%   F   an object of class field representing the scalar field obtained by
%       v.' * voxmfd.g * w
%--------------------------------------------------------------------------
% DEVELOPER TODOs:
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow
%--------------------------------------------------------------------------

%% Main function
%--------------------------------------------------------------------------

F = v.' * ( voxmfd.g * w );
           
return