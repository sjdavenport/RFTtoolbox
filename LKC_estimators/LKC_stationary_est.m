function [ L, L0, nonstatInt ] = LKC_stationary_est( field, dfield, version )
% LKC_stationary_est( cfield, dcfield, version ) estimates the LKCs assuming
% stationarity.
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
%  version a logical/ logical vector. Length depends on voxmfd.D
%          - D = 1, always true.
%          - D = 2, either 1 or 0. 1 if L1 should be estimated, 0 else.
%          - D = 3, logical of length 3. version(1), indicates whether L2
%          should be estimated, version(2) whether the first integral is
%          used in L1 and version(3) whether the second integral is used.
%
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
D = field.D;

if ~exist( 'version', 'var' )
    if field.D < 3
        version = true;
    else
        version = true( [ 1 3 ] );
    end
end

%% Main function
%--------------------------------------------------------------------------

% Construct VoxManifold object by providing Riemannian metric
g = Riemmetric_est( field, dfield );
g = Mask( g );
G = reshape( g.field, [ prod( g.masksize ), D, D ] );
G = constfield( squeeze( mean( G( g.mask(:), :, : ) ) ), g.masksize );
g.field =  G.field;
voxmfd  = VoxManifold( g );

% Obtain the LKCs
[ L, L0, nonstatInt ] = LKC_est( voxmfd, version );

return