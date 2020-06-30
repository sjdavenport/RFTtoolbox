function voxmfd = Field2VoxManifold( field, dfield, d2field )
% ConvField2VoxManifold( field, dfield, masked ) returns the VoxManifold
% object from a convolution field and its derivative. The manifold is given
% by the mask and the Riemannian metric is the induced metric from the
% convolution field.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  field   an object of class Field representing observations of a field,
%          fiberD = 1 and fibersize > 1.
%  dfield  an object of class Field representing observations of the 
%          derivatives of a field.
%
% Optional
%  d2field  an object of class Field representing observations of the 
%           second derivatives of a field.
%
%--------------------------------------------------------------------------
% OUTPUT
% voxmfd  an object of class VoxManifold
%
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% Author: Fabian Telschow
%--------------------------------------------------------------------------

%% Main function
%--------------------------------------------------------------------------

% Construct VoxManifold object by providing Riemannian metric and resadd
% and enlarge
voxmfd = VoxManifold( Riemmetric_est( field, dfield ) );
   
if exist( d2field, 'var' )
    voxmfd.Gamma = Christoffel_est( field, dfield, d2field );
end

return