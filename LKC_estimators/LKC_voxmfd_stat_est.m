function [ L, L0, Lambda ] = LKC_voxmfd_stat_est( field, dfield, version, scale )
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
%                   Default is 1.
%          - D = 3, logical of length 3. version(1), indicates whether L2
%          should be estimated, version(2) whether the first integral is
%          used in L1 and version(3) whether the second integral is used.
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
        version = [true,true,false];
    end
end

if ~exist('scale', 'var')
   scale = 0; 
end


%% Main function
%--------------------------------------------------------------------------
% Construct VoxManifold object by providing Riemannian metric
g = Riemmetric_est( field, dfield );
g = Mask( g );
G = reshape( g.field, [ prod( g.masksize ), D, D ] );

% Calculate the Lambda matrix as the average
Lambda = squeeze( mean( G( g.mask(:), :, : ) ) );

% Adjust by the scaling factor
if scale == 1
    Lambda = Lambda*(field.fibersize-3)/(field.fibersize-2);
end

% Obtain a constant field
G = constfield( Lambda, field.mask );
g.field =  G.field;
voxmfd  = VoxManifold( g );

% Obtain the LKCs
[ L, L0 ] = LKC_est( voxmfd, version );

return
