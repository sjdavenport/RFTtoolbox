function [ L, L0 ] = LKC_est( voxmfd )
% LKC_est( voxmfd ) computes the Lipschitz Killing curvatures for a voxel
% manifold.
%
% While 1D and 2D allow for general non-isotropic fields, the estimate in
% 3D assumes 'local isotropy' for L1. This might be fixed in the future.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  voxmfd  an object of class VoxManifold
%--------------------------------------------------------------------------
% OUTPUT
%   L   a 1 x voxmfd.D vector containing the LKCs L1,..., L_voxmfd.D
%   L0  an integer containing the Euler characteristic of the mask or
%       equivalently the zeroth LKC.
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

% Check that method is implemented for dimension D
if voxmfd.D > 3
    error( strcat( 'D must be < 4. LKC estimation for higher dimensional ',...
                   'voxel manifolds has not been implemented yet.' ) )
end

%% Main function
%--------------------------------------------------------------------------

% Mask the properties of type Field, if not yet done
if ~voxmfd.masked
    voxmfd = Mask( voxmfd );
end

%%% Allocate variables
% Save the voxmfd properties in simpler variables
mask    = voxmfd.mask;
xvals   = voxmfd.xvals;
g       = voxmfd.g.field;
D       = voxmfd.D;
resadd  = voxmfd.resadd;
clear voxmfd

% Allocate vector for Lipschitz Killing curvature
L = NaN * ones( [ 1 D ] );

%%% Get the weights for integration
if resadd ~= 0
    weights = getweights( mask );
    % Reduce weights matrix only to active voxels for speed up
    weights = weights( mask(:) );
else
    weights = ones( [ sum( mask(:) ), 1 ] );
end

%%%%%% BEGIN estimate the LKCs in different dimensions
%%% Compute 0th LKC
L0 = EulerChar( mask, 0.5, D );

%%% Compute LKCs for 0 < d < D
switch D
    case 1
        %%% Calculate LKC1        
        % Voxel size
        xvec = xvals{1};
        
        % Get the volume form
        vol_form = sqrt( max( g(:,1), 0 ) );
        
        % Restrict vol_form and dx to mask
        xvec     = xvec( mask );
        vol_form = vol_form( mask );     

        % Estimate of L1 by integrating volume form over the domain using
        % the trapezoid rule
        L(1) = diff( xvec ) * ( vol_form(1:end-1) + vol_form(2:end) ) / 2;
        
    case 2
        % Get the voxel grid dimensions after resolution increase
        dx = diff( xvals{1} );
        dx = dx(1);
        dy = diff( xvals{2} );
        dy = dy(1);
        
        % Short cuts for the metric entries
        g_xx = g(:, :, 1, 1 );
        g_yy = g(:, :, 2, 2 );
        g_xy = g(:, :, 1, 2 );
        
        %%% Calculate LKC 2
        % Get the volume form, max introduced for stability, since at the
        % boundaries there might be tiny negative numbers due to numerical
        % calculations
        vol_form = sqrt( max( g_xx .* g_yy - g_xy.^2, 0 ) );
        
        % Restrict vol_form to mask
        vol_form = vol_form( mask );
        
        % Integate volume form over the domain. It assumes that each voxel
        % has the same volume dx*dy and simple midpoint integration is used
        % note that we also use weights to give an appropriate quotient to
        % boundary voxels, if resadd > 0, which takes into account that the
        % volume of the boundary voxels needs to be halved, or quartered, etc
        L(2) = sum( vol_form(:) .* weights(:) ) * dx * dy;
            
        %%% Calculate LKC 1
        % Find x and y shift boundary voxels, i.e. horizontal boundary
        % and vertical parts.
        [ bdry, weights] = bndry_voxels( mask, [ "x", "y" ] );
        
        % Integrate using trapozoid rule.
        % Note that we later need to remove half of the end points value,
        % of each line segment which are contained in the x and the y shift
        % boundary. They will be count double otherwise.
        L(1)  = 0.5 *...
                ( sum( sqrt( g_xx( bdry.x ) ) .* weights.x( bdry.x ) ) * dx...
                + sum( sqrt( g_yy( bdry.y ) .* weights.y( bdry.y ) ) ) * dy );
                               
    case 3
        % Get the voxel grid dimensions after resolution increase
        dx = diff( xvals{1} );
        dx = dx(1);
        dy = diff( xvals{2} );
        dy = dy(1);
        dz = diff( xvals{3} );
        dz = dz(1);
        
        % Short cuts for the metric entries
        g_xx = g(:, :, :, 1, 1 );
        g_yy = g(:, :, :, 2, 2 );
        g_zz = g(:, :, :, 3, 3 );
        g_xy = g(:, :, :, 1, 2 );
        g_xz = g(:, :, :, 1, 3 );
        g_yz = g(:, :, :, 2, 3 );
        
        %%% Calculate LKC 3
        % Get the volume form, i.e. sqrt(det g), max introduced for
        % stability
        vol_form = sqrt( max(   g_xx.*g_yy.*g_zz...
                              + g_xy.*g_yz.*g_xz...
                              + g_xz.*g_xy.*g_yz...
                              - g_xz.^2.*g_yy...
                              - g_xy.^2.*g_zz...
                              - g_xx.*g_yz.^2, 0 ) );

        % Restrict vol_form to mask
        vol_form = vol_form( mask );
        
        % Integate volume form over the domain assuming each voxel having
        % the same volume dxdydz. Simple midpoint integration is used.
        L(3) = vol_form(:)' * weights(:) * dx * dy * dz;                  

        %%%%%% Calculate LKC2 and LKC1
        % Find faces and edges and their integration weights. Note that the
        % weights for edges assume that the unit sphere no matter how
        % oriented always has 1/4 or 3/4 of its vectors pointing inside the
        % voxel manifold. I think that should be the case, yet needs to be
        % carefully considered.
        [ bdry, weights ] = bndry_voxels( mask );
        
        %%% Calculate LKC 2;
        % Integrate volume form of faces to get LKC2
        weights.xy = weights.xy( weights.xy ~= 0 );        
        L(2) = sum( sqrt( max( g_xx( bdry.xy ) .* g_yy( bdry.xy )...
                              - g_xy( bdry.xy ).^2, 0 ) )...
                              .* weights.xy(:) ) * dx * dy / 2;
                          
        weights.xz = weights.xz( weights.xz ~= 0 );
        L(2) = L(2) + sum( sqrt( max( g_xx( bdry.xz ) .* g_zz( bdry.xz )...
                              - g_xz( bdry.xz ).^2, 0 ) )...
                              .* weights.xz(:) ) * dx * dz / 2;

        weights.yz = weights.yz( weights.yz ~= 0 );
        L(2) = L(2) + sum( sqrt( max( g_yy( bdry.yz ) .* g_zz( bdry.yz )...
                              - g_yz( bdry.yz ).^2, 0 ) )...
                              .* weights.yz(:) ) * dy * dz / 2;
        
        %%% Calculate LKC 1
        % Integrate volume form of edges to get LKC2
        weights.x = weights.x( weights.x ~= 0 );        
        L(1) = sum( sqrt( max( g_xx( bdry.x ), 0 ) )...
                              .* weights.x(:) ) * dx;
                          
        weights.z = weights.z( weights.z ~= 0 );
        L(1) = L(1) + sum( sqrt( max( g_zz( bdry.z ), 0 ) )...
                              .* weights.z(:) ) * dz;

        weights.y = weights.y( weights.y ~= 0 );
        L(1) = L(1) + sum( sqrt( max( g_yy( bdry.y ), 0 ) )...
                              .* weights.y(:) ) * dy;
end
%%%%%% END estimate the LKCs in different dimensions
           
return