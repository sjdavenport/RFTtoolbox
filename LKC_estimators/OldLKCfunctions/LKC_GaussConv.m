function [L, out] = LKC_GaussConv( Y, nu, D, resAdd, remove )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
%   Y (array T_1 x ... x T_D x N): data array. Last index enumerates the
%                                  samples
%   nu (array 1x1 or 1xD): bandwidth of the Gaussian kernel, if numeric an
%                          isotropic kernel is assumed.  
%   D (array D+1x1): dimension of the domain of the kernel
%   resAdd (integer): denotes the amount of voxels to increase resolution
%   remove (integer): only for theoretical simulations, to account for
%                     boundary effects. Default=0. Don't submit any value,
%                     if you are not simulating theoretical processes.
%
% Output:
%   L (array 1xD): estimated LKCs
%   vol_form (array): values of the estimated volume form from the induced
%                     Gaussian process from Y.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check input and set default values
if nargin < 5
    remove = 0;
end

% Dimension to remove boundary effects
sY     = size( Y ); 
simDim = sY( 1 : D );
N      = sY( D + 1 );

L = NaN * ones( [ 1 D ] );

dx = 1/(resAdd+1);

% get smoothing kernel on higher resolution
[ h, dxh, dyh, dzh ] = gaussNormDerFilt( nu, D, ones([1 D])*dx );

% Dimensions for increased resolution
simDimhr = ( simDim - 1 ) * resAdd + simDim;
   
% number of points which needs to be removed, since they don't belong into
% the estimation regime
remove2 = remove * ( 1 + resAdd);

% sturcture to output the different computed fields for debugging purposes
out = struct();

switch D
    case 1
        % increase the resolution of the raw data by introducing zeros
        Y2 = zeros( [ simDimhr, N ] );
        Y2( 1 : ( resAdd + 1 ) : end, : ) = Y;
        
        % get the convolutional field
        smY = convn(Y2', h, 'same')';
     
        % get the derivative of the convolutional field
        smYx = convn(Y2', dxh, 'same')';
        
        % get the estimates of the covariances
        VY    = var(smY,  0, D+1);
        VdxY  = var(smYx, 0, D+1);
        CYdxY = sum( ( smYx - mean(smYx, D+1) ) .* ...
                     ( smY  - mean(smY,  D+1) ), D+1 ) / (N-1);
                 
        % remove padded values in simulation of generated process to
        % avoid boundary effect
        VY    = VY( (remove2+1):(end-remove2) );
        VdxY  = VdxY( (remove2+1):(end-remove2) );
        CYdxY = CYdxY( (remove2+1):(end-remove2) );
                 
        % get the volume form
        vol_form = sqrt( ( VdxY.*VY - CYdxY.^2 ) ) ./ VY;
        
        % estimate of L1 by integrating volume form over the domain
        L(1) = sum( ( vol_form(1:end-1) + vol_form(2:end) ) ) * dx / 2;
        
        % Fill the output structure
        out.vol_form = vol_form;
        out.VY       = VY;
        out.VdxY     = VdxY;
        out.CYdxY    = CYdxY;

    case 2
        % increase the resolution of the raw data by introducing zeros
        Y2 = zeros( [ simDimhr, N ] );
        Y2( 1:( resAdd + 1 ):end, 1:( resAdd + 1 ):end, : ) = Y;
        
        % get the convolutional field
        smY  = convn( Y2, h, 'same' );
        % get the derivative of the convolutional field
        smYx = convn( Y2, dxh, 'same' );
        smYy = convn( Y2, dyh, 'same' );
        
        % Get the estimates of the covariances
        VY    = var( smY,  0, D+1 );
        VdxY = var( smYx, 0, D+1 );
        VdyY = var( smYy, 0, D+1 );
        CdxYdyY = sum( ( smYy - mean( smYy, D+1 ) ) .* ...
                     ( smYx - mean( smYx, D+1 ) ), D+1 ) / (N-1);
        CYdxY = sum( ( smYx - mean( smYx, D+1 ) ) .* ...
                     ( smY  - mean( smY,  D+1 ) ), D+1 ) / (N-1);
        CYdyY = sum( ( smYy - mean( smYy, D+1 ) ) .* ...
                     ( smY  - mean( smY,  D+1 ) ), D+1 ) / (N-1);
                 
        % entries of riemanian metric
        g_xx = -CYdxY.^2 ./ VY.^2 + VdxY ./ VY;
        g_yy = -CYdyY.^2 ./ VY.^2 + VdyY ./ VY;
        g_xy = -CYdyY .* CYdxY ./ VY.^2 + CdxYdyY ./ VY;
        
        % cut it down to the valid part of the domain
        g_xx = max( g_xx( (remove2 + 1):(end-remove2), (remove2 + 1):(end-remove2) ), 0 );
        g_yy = max( g_yy( (remove2 + 1):(end-remove2), (remove2 + 1):(end-remove2) ), 0 );
        g_xy = g_xy( (remove2 + 1):(end-remove2), (remove2 + 1):(end-remove2) );
        
        % get the volume form, max intorduced for stability
        vol_form = sqrt( max( g_xx.*g_yy - g_xy.*g_xy, 0 ) );

%        figure
%        imagesc(vol_form)
%        colorbar;
        
        %%%% calculate the Lipschitz killing curvatures
        L(1) = sum(...
                    sqrt(g_xx(1,1:end-1)')     + sqrt(g_xx(1,2:end)') + ...
                    sqrt(g_yy(1:end-1,1))      + sqrt(g_yy(2:end,1) ) + ...
                    sqrt(g_xx(end-1,1:end-1)') + sqrt(g_xx(end-1,2:end)' )+ ...
                    sqrt(g_yy(1:end-1,end-1))  + sqrt(g_yy(2:end,end-1) )...
                   ) * dx / 2 / 2;
        
        % get meshgrid of domain and delaunay triangulation for integration
        % over the domain
        [ Xgrid, Ygrid ] = meshgrid( 1:dx:(sY(1)-2*remove),  1:dx:(sY(2)-2*remove) );   
        DT = delaunayTriangulation( [ Xgrid(:), Ygrid(:) ] );
        
        L(2) = integrateTriangulation( DT, vol_form(:) );
        
        % Fill the output structure
        out.vol_form = vol_form;
        out.VY       = VY;
        out.VdxY     = VdxY;
        out.CYdxY    = CYdxY;
        
        out.VdyY    = VdyY;
        out.CYdyY   = CYdyY;
        out.CYdxY   = CYdxY;
        out.CdxYdyY = CdxYdyY;
        
    case 3
        % increase the resolution of the raw data by introducing zeros
        Y2 = zeros( [ simDimhr, N ] );
        Y2( 1:( resAdd + 1 ):end, 1:( resAdd + 1 ):end, : ) = Y;
        
        % get the convolutional field
        smY  = convn( Y2, h, 'same' );
        % get the derivative of the convolutional field
        smYx = convn( Y2, dxh, 'same' );
        smYy = convn( Y2, dyh, 'same' );
        smYz = convn( Y2, dyz, 'same' );
        
        % Get the estimates of the covariances
        VY    = var( smY,  0, D+1 );
        VdxY = var( smYx, 0, D+1 );
        VdyY = var( smYy, 0, D+1 );
        VdzY = var( smYy, 0, D+1 );
        
        CdxYdyY = sum( ( smYy - mean( smYy, D+1 ) ) .* ...
                     ( smYx - mean( smYx, D+1 ) ), D+1 ) / (N-1);
        CdxYdzY = sum( ( smYx - mean( smYx, D+1 ) ) .* ...
                     ( smYz - mean( smYz, D+1 ) ), D+1 ) / (N-1);
        CdyYdzY = sum( ( smYy - mean( smYy, D+1 ) ) .* ...
                     ( smYz - mean( smYz, D+1 ) ), D+1 ) / (N-1);

        CYdxY = sum( ( smYx - mean( smYx, D+1 ) ) .* ...
                     ( smY  - mean( smY,  D+1 ) ), D+1 ) / (N-1);
        CYdyY = sum( ( smYy - mean( smYy, D+1 ) ) .* ...
                     ( smY  - mean( smY,  D+1 ) ), D+1 ) / (N-1);
        CYdzY = sum( ( smYz - mean( smYz, D+1 ) ) .* ...
                     ( smY  - mean( smY,  D+1 ) ), D+1 ) / (N-1);
                 
        % entries of riemanian metric
        g_xx = -CYdxY.^2 ./ VY.^2 + VdxY ./ VY;
        g_yy = -CYdyY.^2 ./ VY.^2 + VdyY ./ VY;
        g_zz = -CYdzY.^2 ./ VY.^2 + VdzY ./ VY;
        g_xy = -CYdyY .* CYdxY ./ VY.^2 + CdxYdyY ./ VY;
        g_xz = -CYdzY .* CYdxY ./ VY.^2 + CdxYdzY ./ VY;
        g_yz = -CYdzY .* CYdyY ./ VY.^2 + CdyYdzY ./ VY;
        
        % cut it down to the valid part of the domain
        g_xx = max( g_xx( (remove2 + 1):(end-remove2), (remove2 + 1):(end-remove2) ), 0 );
        g_yy = max( g_yy( (remove2 + 1):(end-remove2), (remove2 + 1):(end-remove2) ), 0 );
        g_zz = max( g_zz( (remove2 + 1):(end-remove2), (remove2 + 1):(end-remove2) ), 0 );
        g_xy = g_xy( (remove2 + 1):(end-remove2), (remove2 + 1):(end-remove2) );
        g_xz = g_xz( (remove2 + 1):(end-remove2), (remove2 + 1):(end-remove2) );
        g_yz = g_yz( (remove2 + 1):(end-remove2), (remove2 + 1):(end-remove2) );

        % get the volume form, max intorduced for stability
        vol_form = sqrt( max(   g_xx.*g_yy.*g_zz + g_xy.*g_yz.*g_xz + g_xz.*g_xy.*g_yz...
                              - g_xz.*g_yy.*g_xz - g_xy.*g_xy.*g_zz - g_xx.*g_yz.*g_yz, 0 ) );
end    
end