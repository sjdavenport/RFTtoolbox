function voxmfd = ConvField2VoxManifold( cfield, dcfield, masked )
% ConvField2VoxManifold( cfield, dcfield, masked ) returns the VoxManifold
% object from a convolution field and its derivative. The manifold is given
% by the mask and the Riemannian metric is the induced metric from the
% convolution field.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  cfield   an object of class ConvField of derivtype 0, fiberD = 1 and
%           fibersize > 1.
%  dcfield  an object of class ConvField of derivtype 1, fiberD = 1 and
%           fibersize > 1.
% Optional
%  masked   a logical indicating whether the objects of class Field in the
%           VoxManifold output are masked. Default 1, i.e. they are masked.
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

%% Check optional input
%--------------------------------------------------------------------------

if ~exist( 'masked', 'var' )
    masked = 1;
end

%% Main function
%--------------------------------------------------------------------------

% Construct VoxManifold object by providing Riemannian metric and resadd
% and enlarge
voxmfd = VoxManifold( Lambda_est( cfield, dcfield ),...
                    cfield.resadd, cfield.enlarge );
                
% Mask the Field properties in voxmfd, if required.
if masked
    voxmfd = Mask( voxmfd );
end

return