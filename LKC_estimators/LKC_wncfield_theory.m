function Ltheory = LKC_wncfield_theory( FWHM, D, resAdd, mask, pad )
% LKC_GaussConvTheory( FWHM, D, resAdd )
% computes theoretical Lipschitz Killing curvatures for a convolution field
% derived from a seperable kernel with underlying discrete independent
% gaussian white noise process.
% It uses the fact that the voxels are independent and thereby double
% convolutions become single convolutions.
% Currently, only 1D and 2D is implemented. Moreover, the domain of
% the field is considered to be a box.
%
% This function is mainly for comparison purposes.
%
%--------------------------------------------------------------------------
% ARGUMENTS
%   FWHM    array 1x1 or 1xD containing the FWHM for different directions
%           for smoothing with a Gaussian kernel, if numeric an isotropic
%           kernel is assumed
%   D       integer dimension of the domain of the field
%   resAdd  integer denoting the amount of voxels padded between existing
%           voxels to increase resolution
%   mask    boolean array T_1 x ... x T_D indicating which voxels belong to
%           to the region of interest. (currently must be all ones!!!)
%--------------------------------------------------------------------------
% OUTPUT
%   L_bdryEff 1xD array of theoretical LKCs, if the boundary effect is
%   taken into account
%   L_nobdryEff 1xD array of theoretical LKCs without the boundary effect
%--------------------------------------------------------------------------
% EXAMPLES
%--------------------------------------------------------------------------
% AUTHORS: Fabian Telschow
%--------------------------------------------------------------------------

%------------ compute further constants and allocate variables ------------
% dimension of domain
domDim = size( mask );
% stepsize for inbetween voxel resolution increase
dx = 1 / ( resAdd + 1 );

% Dimensions for field with increased resolution
if D == 1
    domDimhr = ( domDim - 1 ) * resAdd + domDim +...
                        ( [ 0 (resAdd + 1)* 2*pad ] );
else
    domDimhr = ( domDim - 1 ) * resAdd + domDim +...
                ( (resAdd + 1)* 2*pad );
end
% number of points which needs to be removed from increased resolution
% image, since they don't belong into the estimation regime
padhr = (resAdd + 1) * pad;

%%%% compute the theoretical LKC value, this assumes that the
%%%% voxels are independent!
switch D
    case 1
        % increase the resolution of the raw data by introducing zeros
        onesField = zeros( domDimhr );
        onesField( 1 : ( resAdd + 1 ) : end ) = 1;
        
        % convolution kernel and derivatives to be used with convn
        if pad==0
           x = -ceil(1.7*FWHM):dx:ceil(1.7*FWHM);
        else
           x = -pad:dx:pad;
        end  
        [ h, dxh ] = Gker( x, FWHM, 1 );   

        VY    = convn( onesField, h.^2, 'same' );
        VdxY  = convn( onesField, dxh.^2, 'same' );
        CYdxY = convn( onesField, dxh.*h, 'same' );

        if pad ~= 0
            VY    = VY( (padhr+1):(end-padhr) );
            VdxY  = VdxY( (padhr+1):(end-padhr) );
            CYdxY = CYdxY( (padhr+1):(end-padhr) );
        end
        
        % get the volume form
        vol_form = sqrt( ( VdxY.*VY - CYdxY.^2 ) ) ./ VY;

%         figure(1)
%         plot(vol_form)

        % estimate of L1 by integrating volume form over the domain
        Ltheory = sum( ( vol_form(1:end-1) + vol_form(2:end) ) ) * dx / 2;

    case 2
        % initialize output array
        Ltheory = NaN * ones( 1, D );
        % increase the resolution of the raw data by introducing zeros
        onesField = zeros( domDimhr );
        onesField( 1:( resAdd + 1 ):end, 1:( resAdd + 1 ):end ) = 1;
        
        % grid for convolution kernel
        if pad == 0
            siz = ceil( 1.7*FWHM );
        else
            siz = pad;
        end
        [x,y] = meshgrid( -siz:dx:siz, -siz:dx:siz );
        xvals = [x(:), y(:)]';
        
        % convolution kernels to be used ith convn
        h   = reshape( GkerMV( xvals, FWHM ), size(x) );
        dh  = GkerMVderiv( xvals, FWHM );
        dxh = reshape( dh(1,:), size(x) );
        dyh = reshape( dh(2,:), size(x) );

        VY      = convn( onesField, h.^2, 'same' );
        VdxY    = convn( onesField, dxh.^2, 'same' );
        VdyY    = convn( onesField, dyh.^2, 'same' );
        CYdyY   = convn( onesField, dyh.*h, 'same' );
        CYdxY   = convn( onesField, dxh.*h, 'same' );
        CdxYdyY = convn( onesField, dxh.*dyh, 'same' );

        % entries of riemanian metric
        g_xx = -CYdxY.^2 ./ VY.^2 + VdxY ./ VY;
        g_yy = -CYdyY.^2 ./ VY.^2 + VdyY ./ VY;
        g_xy = -CYdyY .* CYdxY ./ VY.^2 + CdxYdyY ./ VY;
        
        % cut it down to the valid part of the domain
        if pad ~= 0
            g_xx = max( g_xx( (padhr+1):(end-padhr), (padhr+1):(end-padhr) ), 0 );
            g_yy = max( g_yy( (padhr+1):(end-padhr), (padhr+1):(end-padhr) ), 0 );
            g_xy = g_xy( (padhr+1):(end-padhr), (padhr+1):(end-padhr) );
        end
        % get the volume form, max intorduced for stability
        vol_form = sqrt( max( g_xx.*g_yy - g_xy.*g_xy, 0 ) );

        %%%% calculate the Lipschitz killing curvatures
        Ltheory(1) = sum(...
                    sqrt(g_xx(1,1:end-1)')     + sqrt(g_xx(1,2:end)') + ...
                    sqrt(g_yy(1:end-1,1))      + sqrt(g_yy(2:end,1) ) + ...
                    sqrt(g_xx(end-1,1:end-1)') + sqrt(g_xx(end-1,2:end)' )+ ...
                    sqrt(g_yy(1:end-1,end-1))  + sqrt(g_yy(2:end,end-1) )...
                   ) * dx / 2 / 2;

        % get meshgrid of domain and delaunay triangulation for integration
        % over the domain
        [ Xgrid, Ygrid ] = meshgrid( 1:dx:domDim(1), ...
                                     1:dx:domDim(2) );   
        DT = delaunayTriangulation( [ Xgrid(:), Ygrid(:) ] );

        Ltheory(2) = integrateTriangulation( DT, vol_form(:) );

     case 3
        % initialize output array
        Ltheory = NaN * ones( 1, D );
        % increase the resolution of the raw data by introducing zeros
        onesField = zeros( domDimhr );
        onesField( 1:( resAdd + 1 ):end, 1:( resAdd + 1 ):end ) = 1;
        
        % grid for convolution kernel
        if pad == 0
            siz = ceil( 1.7*FWHM );
        else
            siz = pad;
        end
        [x,y] = meshgrid( -siz:dx:siz, -siz:dx:siz );
        xvals = [x(:), y(:)]';
        
        % convolution kernels to be used ith convn
        h   = reshape( GkerMV( xvals, FWHM ), size(x) );
        dh  = GkerMVderiv( xvals, FWHM );
        dxh = reshape( dh(1,:), size(x) );
        dyh = reshape( dh(2,:), size(x) );

        VY      = convn( onesField, h.^2, 'same' );
        VdxY    = convn( onesField, dxh.^2, 'same' );
        VdyY    = convn( onesField, dyh.^2, 'same' );
        CYdyY   = convn( onesField, dyh.*h, 'same' );
        CYdxY   = convn( onesField, dxh.*h, 'same' );
        CdxYdyY = convn( onesField, dxh.*dyh, 'same' );

        % entries of riemanian metric
        g_xx = -CYdxY.^2 ./ VY.^2 + VdxY ./ VY;
        g_yy = -CYdyY.^2 ./ VY.^2 + VdyY ./ VY;
        g_xy = -CYdyY .* CYdxY ./ VY.^2 + CdxYdyY ./ VY;
        
        % cut it down to the valid part of the domain
        if pad ~= 0
            g_xx = max( g_xx( (padhr+1):(end-padhr), (padhr+1):(end-padhr) ), 0 );
            g_yy = max( g_yy( (padhr+1):(end-padhr), (padhr+1):(end-padhr) ), 0 );
            g_xy = g_xy( (padhr+1):(end-padhr), (padhr+1):(end-padhr) );
        end
        % get the volume form, max intorduced for stability
        vol_form = sqrt( max( g_xx.*g_yy - g_xy.*g_xy, 0 ) );

        %%%% calculate the Lipschitz killing curvatures
        Ltheory(1) = sum(...
                    sqrt(g_xx(1,1:end-1)')     + sqrt(g_xx(1,2:end)') + ...
                    sqrt(g_yy(1:end-1,1))      + sqrt(g_yy(2:end,1) ) + ...
                    sqrt(g_xx(end-1,1:end-1)') + sqrt(g_xx(end-1,2:end)' )+ ...
                    sqrt(g_yy(1:end-1,end-1))  + sqrt(g_yy(2:end,end-1) )...
                   ) * dx / 2 / 2;

        % get meshgrid of domain and delaunay triangulation for integration
        % over the domain
        [ Xgrid, Ygrid ] = meshgrid( 1:dx:domDim(1), ...
                                     1:dx:domDim(2) );   
        DT = delaunayTriangulation( [ Xgrid(:), Ygrid(:) ] );

        Ltheory(2) = integrateTriangulation( DT, vol_form(:) );
end
end