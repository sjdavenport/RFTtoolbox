function voxmfd = ConvField2VoxManifold( cfield, dcfield, masked )
% ConvField2VoxManifold( obj1, obj2 ) constructs a VoxManifold object from
% convolution fields.
% The convolution fields should be generated using convfield.m.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  cfield   an object of class ConvField of derivtype 0 and fiberD = 1 and
%           fibersize > 1.
%  dcfield  an object of class ConvField of derivtype 1 and fiberD = 1 and
%           fibersize > 1.
% Optional
%  masked   a logical indicating whether the objects of class Field in the
%           VoxManifold output are masked. Default 1, i.e. they are masked.
%
%--------------------------------------------------------------------------
% OUTPUT
% obj  an object of class VoxManifold
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