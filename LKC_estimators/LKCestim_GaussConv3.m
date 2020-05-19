function [L, geom] = LKCestim_GaussConv3( Y, FWHM, mask, resAdd, remove )
% LKCestim_GaussConv( Y, nu, mask, resAdd, remove )
% estimates the Lipschitz Killing curvatures for a convolution process.
% It uses the fact that derivatives can be represented as convolutions
% with the derivative kernel.
% Currently, only 1D and 2D are tested and working. Moreover, the domain of
% the field is considered to be a box.
%
% Required additions:
%   - allow voxel dimensions to be different in different directions,
%   currently distance between voxels is assumed to be 1
%   - include mask
%   - add full 3D estimation
%
%--------------------------------------------------------------------------
% ARGUMENTS
%   Y       data array T_1 x ... x T_D x N. Last index enumerates the
%           samples. Note that N > 1 is required!
%   FWHM    array 1x1 or 1xD containing the FWHM for different directions
%           for smoothing with a Gaussian kernel, if numeric an isotropic
%           kernel is assumed
%   mask    boolean array T_1 x ... x T_D indicating which voxels belong to
%           to the region of interest. (not yet implemented!)
%   resAdd  integer denoting the amount of voxels padded between existing
%           voxels to increase resolution
%   remove  (only for theoretical simulations) integer amount of boundary
%           voxels, which are removed from the boundary.
%           This is only neccessary to handle boundary effects in
%           simulations. Default = 0 (no voxels removed).
%           Only touch, if you are simulating theoretical processes and you
%           want to compare to processes derived from an enlarged domain.
%--------------------------------------------------------------------------
% OUTPUT
%   L       1xD array of estimated LKCs
%   geom    structure containing geometric quantities of the induced metric
%           of the random field.
%--------------------------------------------------------------------------
% EXAMPLES
% %1D
% rf   = noisegen( [35 35], 50, 6 );
% mask = ones([35 35);
% L = LKCestim_GaussConv( rf, 3, mask, 1 );
%
% %2D
% rf   = noisegen( [35 35], 50, 6 );
% mask = ones([35 35);
% L = LKCestim_GaussConv( rf, 3, mask, 1 );
% 
% %3D
% thresh = 1;
% sims = randn(10,10,10)
% clusters = clusterloc(sims, thresh);
% sims > thresh
% clusters
%--------------------------------------------------------------------------
% AUTHORS: Fabian Telschow
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------ get parameters from the random field input ------------------
% Dimension of the input
sY     = size( Y );
% Dimension of the domain
domDim = sY( 1 : end-1 );
% dimension of the domain
D = length( domDim );
% number of samples
nsubj = sY( end );

%------------ check input and set default values --------------------------
if nargin < 5
    remove = 0;
end

% Check mask input (this need to be coded carefully until know only boxes are )
mask = logical(ones(domDim));


%------------ compute further constants and allocate variables ------------
% stepsize for inbetween voxel resolution increase
dx = 1/(resAdd+1);
% Dimensions for field with increased resolution
domDimhr = ( domDim - 1 ) * resAdd + domDim;
% number of points which needs to be removed from increased resolution
% image, since they don't belong into the estimation regime
remove2 = remove * ( 1 + resAdd);

% allocate vector for Lipschitz Killing curvature
L = NaN * ones( [ 1 length(domDim) ] );
% sturcture to output the different computed fields for debugging purposes
geom = struct();

% range for which the kernel is almost zero. kernel will be truncated here.
siz = ceil( 1.7 * FWHM );

%------------ estimate the LKCs -------------------------------------------
switch D
    case 1
        % increase the resolution of the raw data by introducing zeros
        Y2 = zeros( [ domDimhr, nsubj ] );
        Y2( 1:(resAdd + 1):end, : ) = Y;
        
        % grid for convolution kernel
        x   = -siz:dx:siz;
        
        % convolution kernel and derivatives to be used with convn
        [ h, dxh ] = Gker( x, FWHM, 1 );     
        
        % get the convolutional field
        smY = convn( Y2', h, 'same' )';
     
        % get the derivative of the convolutional field
        smYx = convn( Y2', dxh, 'same' )';
        
        % get the estimates of the covariances
        VY    = var( smY,  0, D+1 );
        VdxY  = var( smYx, 0, D+1 );
        CYdxY = sum( ( smYx - mean(smYx, D+1) ) .* smY, D+1 ) / (nsubj-1);
                 
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
        geom.vol_form = vol_form;
        geom.VY       = VY;
        geom.VdxY     = VdxY;
        geom.CYdxY    = CYdxY;

    case 2
        % increase the resolution of the raw data by introducing zeros
        Y2 = zeros( [ domDimhr, nsubj ] );
        Y2( 1:( resAdd + 1 ):end, 1:( resAdd + 1 ):end, : ) = Y;
        
        % grid for convolution kernel
        [x,y] = meshgrid( -siz:dx:siz, -siz:dx:siz );
        xvals = [x(:), y(:)]';
        
        % convolution kernels to be used ith convn
        h   = reshape( GkerMV( xvals, FWHM ), size(x) );
        dh  = GkerMVderiv( xvals, FWHM );
        dxh = reshape( dh(1,:), size(x) );
        dyh = reshape( dh(2,:), size(x) );
        
        % get the convolutional field
        smY  = convn( Y2(:,:,:,1), h, 'same' );
        % get the derivative of the convolutional field
        smYx = convn( Y2, dxh, 'same' );
        smYy = convn( Y2, dyh, 'same' );
        
        % Get the estimates of the covariances
        VY   = var( smY,  0, D+1 );
        VdxY = var( smYx, 0, D+1 );
        VdyY = var( smYy, 0, D+1 );
        CdxYdyY = sum( ( smYy - mean( smYy, D+1 ) ) .* ...
                         smYx, D+1 ) / (nsubj-1);
        CYdxY = sum( ( smYx - mean( smYx, D+1 ) ) .* ...
                       smY, D+1 ) / (nsubj-1);
        CYdyY = sum( ( smYy - mean( smYy, D+1 ) ) .* ...
                       smY, D+1 ) / (nsubj-1);
                 
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
        
        %%%% calculate the Lipschitz killing curvatures
        L(1) = (...
                   sum( sqrt(g_xx(1,1:end-1)')     + sqrt(g_xx(1,2:end)') ) + ...
                   sum( sqrt(g_yy(1:end-1,1))      + sqrt(g_yy(2:end,1) ) ) + ...
                   sum( sqrt(g_xx(end-1,1:end-1)') + sqrt(g_xx(end-1,2:end)' ) )+ ...
                   sum( sqrt(g_yy(1:end-1,end-1))  + sqrt(g_yy(2:end,end-1) ) )...
                   ) * dx / 2 / 2;
        
        % get meshgrid of domain and delaunay triangulation for integration
        % over the domain
        [ Xgrid, Ygrid ] = meshgrid( 1:dx:(sY(1)-2*remove), ...
                                     1:dx:(sY(2)-2*remove) );   
        DT = delaunayTriangulation( [ Xgrid(:), Ygrid(:) ] );
        
        L(2) = integrateTriangulation( DT, vol_form(:) );
        
        % Fill the output structure
        geom.vol_form = vol_form;
        geom.VY       = VY;
        geom.VdxY     = VdxY;
        geom.CYdxY    = CYdxY;
        
        geom.VdyY    = VdyY;
        geom.CYdyY   = CYdyY;
        geom.CYdxY   = CYdxY;
        geom.CdxYdyY = CdxYdyY;
        
    case 3
        % increase the resolution of the raw data by introducing zeros
        Y2 = zeros( [ domDimhr, nsubj ] );
        Y2( 1:( resAdd + 1 ):end, 1:( resAdd + 1 ):end,...
            1:( resAdd + 1 ):end, : ) = Y;
        
        % grid for convolution kernel
        [x,y,z] = meshgrid( -siz:dx:siz, -siz:dx:siz, -siz:dx:siz );
        xvals = [x(:), y(:), z(:)]';
        
        % convolution kernels to be used with convn
        h   = reshape( GkerMV( xvals, FWHM ), size(x) );
        dh  = GkerMVderiv( xvals, FWHM );
        dxh = reshape( dh(1,:), size(x) );
        dyh = reshape( dh(2,:), size(x) );
        dzh = reshape( dh(3,:), size(x) );
        
        % get the convolutional field
        smY  = convn( Y2, h, 'same' );
        % get the derivative of the convolutional field
        smYx = convn( Y2, dxh, 'same' );
        smYy = convn( Y2, dyh, 'same' );
        smYz = convn( Y2, dzh, 'same' );
        
        % Get the estimates of the covariances
        VY   = var( smY,  0, D+1 );
        VdxY = var( smYx, 0, D+1 );
        VdyY = var( smYy, 0, D+1 );
        VdzY = var( smYy, 0, D+1 );
        
        CdxYdyY = sum( ( smYy - mean( smYy, D+1 ) ) .* ...
                     ( smYx - mean( smYx, D+1 ) ), D+1 ) / (nsubj-1);
        CdxYdzY = sum( ( smYx - mean( smYx, D+1 ) ) .* ...
                     ( smYz - mean( smYz, D+1 ) ), D+1 ) / (nsubj-1);
        CdyYdzY = sum( ( smYy - mean( smYy, D+1 ) ) .* ...
                     ( smYz - mean( smYz, D+1 ) ), D+1 ) / (nsubj-1);

        CYdxY = sum( ( smYx - mean( smYx, D+1 ) ) .* ...
                     ( smY  - mean( smY,  D+1 ) ), D+1 ) / (nsubj-1);
        CYdyY = sum( ( smYy - mean( smYy, D+1 ) ) .* ...
                     ( smY  - mean( smY,  D+1 ) ), D+1 ) / (nsubj-1);
        CYdzY = sum( ( smYz - mean( smYz, D+1 ) ) .* ...
                     ( smY  - mean( smY,  D+1 ) ), D+1 ) / (nsubj-1);
                 
        % entries of riemanian metric/ Lambda matrix from neuroimaging
        g_xx = ( -CYdxY.^2 + VdxY .* VY ) ./ VY.^2;
        g_yy = ( -CYdyY.^2 + VdyY .* VY ) ./ VY.^2;
        g_zz = ( -CYdzY.^2 + VdzY .* VY ) ./ VY.^2;
        g_xy = ( -CYdyY .* CYdxY + CdxYdyY .* VY ) ./ VY.^2;
        g_xz = ( -CYdzY .* CYdxY + CdxYdzY .* VY ) ./ VY.^2;
        g_yz = ( -CYdzY .* CYdyY + CdyYdzY .* VY ) ./ VY.^2;
        
        % cut it down to the valid part of the domain
        g_xx = max( g_xx( (remove2 + 1):(end-remove2),...
                          (remove2 + 1):(end-remove2),...
                          (remove2 + 1):(end-remove2) ), 0 );
        g_yy = max( g_yy( (remove2 + 1):(end-remove2),...
                          (remove2 + 1):(end-remove2),...
                          (remove2 + 1):(end-remove2) ), 0 );
        g_zz = max( g_zz( (remove2 + 1):(end-remove2),...
                          (remove2 + 1):(end-remove2),...
                          (remove2 + 1):(end-remove2) ), 0 );
        g_xy = g_xy( (remove2 + 1):(end-remove2),...
                     (remove2 + 1):(end-remove2),...
                     (remove2 + 1):(end-remove2) );
        g_xz = g_xz( (remove2 + 1):(end-remove2),...
                     (remove2 + 1):(end-remove2),...
                     (remove2 + 1):(end-remove2) );
        g_yz = g_yz( (remove2 + 1):(end-remove2),...
                     (remove2 + 1):(end-remove2),...
                     (remove2 + 1):(end-remove2) );
        %------------------------------------------------------------------
        %   compute L3
        %------------------------------------------------------------------
        % get the volume form, max intorduced for stability
        vol_form = sqrt( max(   g_xx.*g_yy.*g_zz...
                              + g_xy.*g_yz.*g_xz...
                              + g_xz.*g_xy.*g_yz...
                              - g_xz.^2.*g_yy...
                              - g_xy.^2.*g_zz...
                              - g_xx.*g_yz.^2, 0 ) );
                          
        % get meshgrid of domain and delaunay triangulation for integration
        % over the domain
        [ Xgrid, Ygrid, Zgrid ] = meshgrid( 1:dx:(sY(1)-2*remove), ...
                                     1:dx:(sY(2)-2*remove), ...
                                     1:dx:(sY(3)-2*remove) );   
        DT = delaunayTriangulation( [ Xgrid(:), Ygrid(:), Zgrid(:) ] );
        
        L(3) = integrateTriangulation( DT, vol_form(:) );
        %------------------------------------------------------------------
        %   compute L2
        %------------------------------------------------------------------
        sG = size(g_xx);
        % compute the volume form of the faces
        ind_xy_b = { ':', ':', 1 };
        ind_xy_t = { ':', ':', sG(3) };
        ind_xz_b = { ':', 1, ':' };
        ind_xz_t = { ':', sG(3), ':' };
        ind_yz_b = { 1 ,':', ':' };
        ind_yz_t = { sG(3), ':', ':' };
        
        % integrate xy_b face volume form
        vol_form = sqrt( max( g_xx(ind_xy_b{:}).*g_yy(ind_xy_b{:})...
                              - g_xy(ind_xy_b{:}).^2, 0 ) );
        [ Xgrid, Ygrid ] = meshgrid( 1:dx:(sY(1)-2*remove), ...
                                     1:dx:(sY(2)-2*remove) ); 
        DT = delaunayTriangulation( [ Xgrid(:), Ygrid(:) ] );
        L(2) = integrateTriangulation( DT, vol_form(:) );
        
        % integrate xy_t face volume form
        vol_form = sqrt( max( g_xx(ind_xy_t{:}).*g_yy(ind_xy_t{:})...
                              - g_xy(ind_xy_t{:}).^2, 0 ) );
        L(2) = L(2) + integrateTriangulation( DT, vol_form(:) );
        
        % integrate xz_b face volume form
        vol_form = sqrt( max( g_xx(ind_xz_b{:}).*g_zz(ind_xz_b{:})...
                              - g_xz(ind_xz_b{:}).^2, 0 ) );
        [ Xgrid, Zgrid ] = meshgrid( 1:dx:(sY(1)-2*remove), ...
                                     1:dx:(sY(3)-2*remove) ); 
        DT = delaunayTriangulation( [ Xgrid(:), Zgrid(:) ] );
        L(2) = L(2) + integrateTriangulation( DT, vol_form(:) );
        
        % integrate xz_t face volume form
        vol_form = sqrt( max( g_xx(ind_xz_t{:}).*g_zz(ind_xz_t{:})...
                              - g_xz(ind_xz_t{:}).^2, 0 ) );
        L(2) = L(2) + integrateTriangulation( DT, vol_form(:) );
        
        % integrate yz_b face volume form
        vol_form = sqrt( max( g_yy(ind_yz_b{:}).*g_zz(ind_yz_b{:})...
                              - g_yz(ind_yz_b{:}).^2, 0 ) );
        [ Ygrid, Zgrid ] = meshgrid( 1:dx:(sY(2)-2*remove), ...
                                     1:dx:(sY(3)-2*remove) ); 
        DT = delaunayTriangulation( [ Ygrid(:), Zgrid(:) ] );
        L(2) = L(2) + integrateTriangulation( DT, vol_form(:) );
        
        % integrate yz_t face volume form
        vol_form = sqrt( max( g_yy(ind_yz_t{:}).*g_zz(ind_yz_t{:})...
                              - g_yz(ind_yz_t{:}).^2, 0 ) );
        L(2) = L(2) + integrateTriangulation( DT, vol_form(:) );
        
        %------------------------------------------------------------------
        %   compute L1 by regression from GKF and plugging in the already
        %   estimated L2 and L3
        %------------------------------------------------------------------
        % probably not necessary, since it might not perform better than
        % the simple HPE.
end    
end