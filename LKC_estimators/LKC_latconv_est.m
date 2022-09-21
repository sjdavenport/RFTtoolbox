function [ L, L0, nonstatInt, cfields ] = LKC_latconv_est( lat_data,...
                                                            params, version )
% LKC_latconv_est( lat_data, params, version ) estimates the LKCs of the
% voxel manifold defined by the domain (aka its mask) of a convolution
% field and its induced Riemannian metric.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  lat_data  an object of class Field representing lattice observations of
%            a field, fiberD = 1 and fibersize > 1.
%  params    an object of class ConvFieldParams.      
%
% Optional
%  version a logical/ logical vector. Length depends on voxmfd.D
%          - D = 1, always true.
%          - D = 2, either 1 or 0. 1 if L1 should be estimated, 0 else.
%          - D = 3, logical of length 3. version(1) indicates whether L2
%                   should be estimated, version(2) whether the first 
%                   integral is used in L1 and version(3) whether also
%                   the second integral is used in L1. Default: [1 1 0];
%                   i.e., the stationary approximation of L1
%
%--------------------------------------------------------------------------
% OUTPUT
% L   an 1 x field.D vector containing the LKCs L_1,...,L_field.D
% L0  a numeric containing the zeroth LKC, i.e. the Euler characteristic of
%     the domain.
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% Author: Fabian Telschow
%--------------------------------------------------------------------------

%% Check and add optional input
%--------------------------------------------------------------------------

if ~exist( 'version', 'var' )
    if lat_data.D < 3
        version = true;
    else
        version = logical( [ 1 1 0 ] );
    end
end

%% Main function
%--------------------------------------------------------------------------
% Compute the convolution fields
cfields    = cell( [ 1 3 ] );
[ tmp, ss ] = convfield( lat_data, params, 0 );
cfields{1} = tmp ./ ss;
cfields{2} = convfield( lat_data, params, 1 )./ ss;
cfields{3} = Field();

% Compute the LKCs
if lat_data.D == 3
    if version(3) == 1
        cfields{3} = convfield( lat_data, params, 2 )./ ss;
        [ L, L0, nonstatInt ]  = LKC_voxmfd_est( cfields{1}, cfields{2},...
                                                 cfields{3}, version );
    else
        [ L, L0, nonstatInt ] = LKC_voxmfd_est( cfields{1}, cfields{2},...
                                                cfields{3}, version );
    end
else
	[ L, L0, nonstatInt ] = LKC_voxmfd_est( cfields{1}, cfields{2},...
                                            cfields{3}, version );
end

return