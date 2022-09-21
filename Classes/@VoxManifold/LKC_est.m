function [ L, L0, nonstatInt ] = LKC_est( voxmfd, version )
% LKC_est( voxmfd, version ) computes the Lipschitz Killing curvatures for 
% a voxel manifold for an arbitrary metric g.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  voxmfd  an object of class VoxManifold. voxmfd.D < 4.
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
%--------------------------------------------------------------------------
% OUTPUT
%   L   a 1 x voxmfd.D vector containing the LKCs L1,..., L_voxmfd.D
%   L0  an integer containing the Euler characteristic of the mask or
%       equivalently the zeroth LKC.
%--------------------------------------------------------------------------
% DEVELOPER TODOs:
%--------------------------------------------------------------------------
% EXAMPLES
% See test_LKC_est.m
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

%% Check optional input and get important constants
%--------------------------------------------------------------------------
if ~exist( 'version', 'var' )
    if voxmfd.D < 3
        version = true;
    else
        version = logical( [ 1 1 0 ] );
    end
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

% Allocate vector for Lipschitz Killing curvature
L = NaN * ones( [ 1 D ] );

%%% Get the weights for integration
weights = getweights( mask );
% Reduce weights matrix only to active voxels for speed up
weights = weights( mask(:) );

%%%%%% BEGIN estimate the LKCs in different dimensions
%%% Compute 0th LKC
L0 = EulerChar( mask, 0.5, D );

%%% Preallocate output for the nonstationary part of L1
nonstatInt = NaN;

%%% Compute LKCs for 0 < d < D
switch D
    case 1
        %%% Calculate LKC1        
        % Voxel size
        dx = diff( xvals{1} );
        dx = dx(1);
        
        % Get the volume form
        vol_form = sqrt( max( g(:,1), 0 ) );

        % Estimate of L1 by integrating volume form over the domain using
        % the trapezoid rule
        L(1) = sum( dx * weights * vol_form( mask ) );
        
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
        if version(1)
            % Find x and y shift boundary voxels, i.e. horizontal boundary
            % and vertical parts.
            [ bdry, weights] = bndry_voxels( mask, [ "x", "y" ] );

            % Integrate using trapozoid rule.
            % Note that we later need to remove half of the end points value,
            % of each line segment which are contained in the x and the y shift
            % boundary. They will be count double otherwise.
            L(1)  = 0.5 *...
                    ( sum( sqrt( g_xx( bdry.x ) ) .* weights.x( bdry.x ) ) * dx...
                    + sum( sqrt( g_yy( bdry.y ) ) .* weights.y( bdry.y ) ) * dy );
        end
                               
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
        vol_form = sqrt( max(   g_xx .* g_yy .* g_zz...
                              + g_xy .* g_yz .* g_xz...
                              + g_xz .* g_xy .* g_yz...
                              - g_xz.^2 .* g_yy...
                              - g_xy.^2 .* g_zz...
                              - g_xx .* g_yz.^2, 0 ) );

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
        [ bdry, weights, eucangle, orientation ] = bndry_voxels( mask );
        
        %%% Calculate LKC 2;
        % Weights for integration
        weights.xy = weights.xy( weights.xy ~= 0 );
        weights.xz = weights.xz( weights.xz ~= 0 );
        weights.yz = weights.yz( weights.yz ~= 0 );
        
        % Volume form on boundaries
        vol_xy = sqrt( max( g_xx( bdry.xy ) .* g_yy( bdry.xy )...
                              - g_xy( bdry.xy ).^2, 0 ) );
        vol_xz = sqrt( max( g_xx( bdry.xz ) .* g_zz( bdry.xz )...
                              - g_xz( bdry.xz ).^2, 0 ) );
        vol_yz = sqrt( max( g_yy( bdry.yz ) .* g_zz( bdry.yz )...
                              - g_yz( bdry.yz ).^2, 0 ) );
        if version(1)                  
            % Integrate volume form of faces to get LKC2
            L(2) =   sum( vol_xy .* weights.xy(:) ) * dx * dy / 2 ...                        
                   + sum( vol_xz .* weights.xz(:) ) * dx * dz / 2 ...
                   + sum( vol_yz .* weights.yz(:) ) * dy * dz / 2;
        end
        
        %%% Calculate LKC 1
        % Get the opening angles with respect to the metric g
        angle = IntegralAlpha( voxmfd, bdry, eucangle );
        angle.x = eucangle.x( eucangle.x ~= 0 );
        angle.x( angle.x > pi ) = - pi/2;
        angle.y = eucangle.y( eucangle.y ~= 0 );
        angle.y( angle.y > pi ) =  - pi/2;      
        angle.z = eucangle.z( eucangle.z ~= 0 );
        angle.z( angle.z > pi ) = - pi/2;

        clear eucangle

        % Integrate volume form of edges against the internal cutting
        % angles with orthogonal plane
        if version(2)
            weights.x = weights.x( weights.x ~= 0 );
            weights.y = weights.y( weights.y ~= 0 );
            weights.z = weights.z( weights.z ~= 0 );
            L(1) = sum( sqrt( max( g_xx( bdry.x ), 0 ) )...
                                  .* weights.x .* angle.x ) * dx;
            L(1) = L(1) + sum( sqrt( max( g_yy( bdry.y ), 0 ) )...
                                  .* weights.y .* angle.y ) * dy;       
            L(1) = L(1) + sum( sqrt( max( g_zz( bdry.z ), 0 ) )...
                                  .* weights.z .* angle.z ) * dz;
        end
                       
        % Integrate volume form of faces versus the trace of the shape operator
        if version(3)
            if isnan( L(1) )
                L(1) = 0;
            end
            
            count = 0;
            
            for type = ["xy", "xz", "yz"]
                count = count + 1;
                % Order needs to ensure that E_k x E_l = E_m
                switch type
                    case "xy"
                        k = 1;
                        l = 2;
                        direc = [ k l 3];
                        dVol = vol_xy .* dx * dy;
                    case "xz"
                        k = 3;
                        l = 1;
                        direc = [ k l 2];
                        dVol = vol_xz .*dx * dz;
                    case "yz"
                        k = 2;
                        l = 3;
                        direc = [ k l 1];
                        dVol = vol_yz .* dy * dz;
                end
                
                % Get the orthonormal frame with first two vectors embedded
                % into the type-face and picking the relevant coordinates
                [ U, V, W ] = OrthNormFrame( voxmfd, direc, bdry.(type), 1 );
                Uk = U( :, k );
                clear U
                Vk = V( :, k );
                Vl = V( :, l );
                clear V
                % Make W outward pointing normal
                W = W .* orientation.(type)( bdry.(type)(:) );

                % Get the vectors derived from the Christoffel symbols
                Gammakk = squeeze( collapse( ...
                                voxmfd.Gamma( :, :, :, k, k, : )...
                           ) );
                Gammakk = Gammakk( bdry.(type)(:), : );
                Gammall = squeeze( collapse( ...
                                voxmfd.Gamma( :, :, :, l, l, : )...
                           ) );
                Gammall = Gammall( bdry.(type)(:), : );
                Gammakl = squeeze( collapse( ...
                                voxmfd.Gamma( :, :, :, k, l, : )...
                           ) );
                Gammakl = Gammakl( bdry.(type)(:), : );

                % Add the integral part belonging to type-face integrals
                integrand = ( ( Uk.^2 + Vk.^2 ) .* ( W.' *  Gammakk ) ...
                                + Vl.^2 .* ( W.' * Gammall ) ...
                                + Vk .* Vl .* ( W.' * Gammakl ) )...
                                .* dVol;
                nonstatInt = sum( integrand.field );
                L(1)       = L(1) +  nonstatInt;
            end
        end
        
        % Get the correct scaling
        L(1) = L(1) / 2 / pi;

end
%%%%%% END estimate the LKCs in different dimensions
           
return