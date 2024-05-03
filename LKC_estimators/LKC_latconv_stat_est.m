function [ L, L0, cfields ] = LKC_latconv_stat_est( lat_data, params, scale )
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
        version = true( [ 1 3 ] );
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
[ L, L0, ~ ] = LKC_voxmfd_stat_est( cfields{1}, cfields{2}, version, scale );

return