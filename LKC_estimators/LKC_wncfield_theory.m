function L = LKC_wncfield_theory( mask, Kernel, resAdd, mask_opt, enlarge )
% computes theoretical Lipschitz Killing curvatures for a convolution field
% derived from a seperable kernel with underlying discrete independent
% Gaussian white noise process.
% It uses the fact that the voxels are independent and thereby double
% convolutions become single convolutions.
%
% This function is mainly for comparison purposes.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   mask      a logical array of dimension T_1 x...x T_D.
%   Kernel    array 1x1 or 1xD containing the FWHM for different directions
%             for smoothing with a Gaussian kernel, if numeric an isotropic
%             kernel is assumed.
% Optional
%   resAdd     integer denoting the amount of voxels padded between 
%              existing voxels to increase resolution. Default 1.
%   mask_opt   2 x 1 logical vector. FIRST COMPONENT, if "1" it applies
%              the mask prior of application of convolution fields.
%              SECOND COMPONENT, if "1" mask is applied after
%              computing geometric properties.
%              Note that [1 1] is possible which means mask is applied in
%              both stages.
%   enlarge    an integer denoting the amount of voxels the resolution
%              increased mask is enlarged by dilation. Note that for
%              unbiased estimation resAdd needs to be an odd number and
%              enlarge needs to be set to the default value currently.
%              Default ceil( resAdd / 2 )
%--------------------------------------------------------------------------
% OUTPUT
%   L  1 x D array of theoretical LKCs computed assuming independence
%      of voxels mean zero and variance 1
%--------------------------------------------------------------------------
% EXAMPLES
%--------------------------------------------------------------------------
% AUTHORS: Fabian Telschow
%--------------------------------------------------------------------------

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Check input and get important constants from the mandatory input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get size of the mask
sM = size( mask );

% Design question, do we want to force the user to use logicals?
mask = logical( mask );

% dimension of the domain, since matlab is not consistent for D<1, we need
% to catch this case 
if sM( 2 ) == 1 && length( sM ) == 2
    D = 1;
else
    D = length( sM );
end

% check that method is implemented for dimension D
if D > 3
    error( strcat( 'D must be < 4. Higher dimensional domains have',...
                   'not been implemented' ) );
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% add/check optional values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist( 'resAdd', 'var' )
   % default number of resolution increasing voxels between observed voxels
   resAdd = 1;
end

if ~exist( 'enlarge', 'var' )
   % default method for enlarge variable
   enlarge = ceil( resAdd / 2 );
end

if ~exist( 'mask_opt', 'var' )
   % default method for mask_opt, which controls when the mask is applied
   % prior/after application of smoothing using convfield.m
   mask_opt = [ 1 1 ];
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  main function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% allocate vector for Lipschitz Killing curvature
L = NaN * ones( [ 1 D ] );

% get mask on higher resolution and the weights of each voxel for the
% volume computation.
if resAdd ~= 0
    [ mask, weights ] = mask_highres( mask, resAdd, enlarge );
    % reduce weights matrix only to active voxels for speed up
    weights = weights( mask );
else
    weights = ones( [ sum( mask(:) ), 1 ] );
end

% get the size of the resolution increased domain
Dimhr = size( mask );


%%%%%% BEGIN get the Riemannian metric/Lambda matrix assuming that the
%%%%%% voxels are independent Gaussian white noise
% preallocate the Riemannian metric
g = NaN * ones( [ Dimhr D D ] );

% create index to fill the original mask at the correct voxels of the high
% resolution mask
index = cell( [ 1 D ] );
for d = 1:D
    index{d} = ( enlarge + 1 ):( resAdd + 1 ):( Dimhr(d) - enlarge );
end

% fill the resolution increased mask with the values from the original mask
% at the correct locations
onesField = zeros( Dimhr );
onesField( index{:} ) = 1;

%%% get the entries of the Riemannian metric
switch D
    case 1
        % get the theoretical variance of the field and the variance of
        % derivatives
        VY    = convn( onesField, h.^2, 'same' );
        VdxY  = convn( onesField, dxh.^2, 'same' );
        CYdxY = convn( onesField, dxh.*h, 'same' );
        
        % get the volume form
        g = ( VdxY.*VY - CYdxY.^2 ) ./ VY;

    case 2
        % get the theoretical variance of the field and the variance of
        % derivatives
        VY      = convn( onesField, h.^2, 'same' );
        VdxY    = convn( onesField, dxh.^2, 'same' );
        VdyY    = convn( onesField, dyh.^2, 'same' );
        CYdyY   = convn( onesField, dyh.*h, 'same' );
        CYdxY   = convn( onesField, dxh.*h, 'same' );
        CdxYdyY = convn( onesField, dxh.*dyh, 'same' );

        % entries of riemanian metric
        g( :, :, 1, 1 ) = -CYdxY.^2 ./ VY.^2 + VdxY ./ VY;
        g( :, :, 2, 2 ) = -CYdyY.^2 ./ VY.^2 + VdyY ./ VY;
        g( :, :, 1, 2 ) = -CYdyY .* CYdxY ./ VY.^2 + CdxYdyY ./ VY;

     case 3
        % Get the estimates of the covariances
        VY   = var( convY,  0, D+1 );
        VdxY = var( convYx, 0, D+1 );
        VdyY = var( convYy, 0, D+1 );
        VdzY = var( convYy, 0, D+1 );
        
        CdxYdyY = sum( ( convYy - mean( convYy, D+1 ) ) .* ...
                       ( convYx - mean( convYx, D+1 ) ), D+1 ) / (nsubj-1);
        CdxYdzY = sum( ( convYx - mean( convYx, D+1 ) ) .* ...
                       ( convYz - mean( convYz, D+1 ) ), D+1 ) / (nsubj-1);
        CdyYdzY = sum( ( convYy - mean( convYy, D+1 ) ) .* ...
                       ( convYz - mean( convYz, D+1 ) ), D+1 ) / (nsubj-1);

        CYdxY = sum( ( convYx - mean( convYx, D+1 ) ) .* ...
                     ( convY  - mean( convY,  D+1 ) ), D+1 ) / (nsubj-1);
        CYdyY = sum( ( convYy - mean( convYy, D+1 ) ) .* ...
                     ( convY  - mean( convY,  D+1 ) ), D+1 ) / (nsubj-1);
        CYdzY = sum( ( convYz - mean( convYz, D+1 ) ) .* ...
                     ( convY  - mean( convY,  D+1 ) ), D+1 ) / (nsubj-1);
                 
        % entries of riemanian metric/ Lambda matrix from neuroimaging
        g( :, :, :, 1, 1 ) = ( -CYdxY.^2 + VdxY .* VY ) ./ VY.^2;
        g( :, :, :, 2, 2 ) = ( -CYdyY.^2 + VdyY .* VY ) ./ VY.^2;
        g( :, :, :, 3, 3 ) = ( -CYdzY.^2 + VdzY .* VY ) ./ VY.^2;
        g( :, :, :, 1, 2 ) = ( -CYdyY .* CYdxY + CdxYdyY .* VY ) ./ VY.^2;
        g( :, :, :, 1, 3 ) = ( -CYdzY .* CYdxY + CdxYdzY .* VY ) ./ VY.^2;
        g( :, :, :, 2, 3 ) = ( -CYdzY .* CYdyY + CdyYdzY .* VY ) ./ VY.^2;
        
end


%%%%%% BEGIN estimate the LKCs in different dimensions
%%% Compute 0th LKC
L0 = EulerChar( mask, 0.5, D );

%%% compute LKCs for 0 < d <= D
switch D
    case 1
        %%% calculate LKC1        
        % voxel size
        xvec = xvals{1};
        
        % get the volume form
        vol_form = sqrt( max( g(:,1), 0 ) );
        
        % restrict vol_form and dx to mask
        if mask_opt(2) == 1
            xvec     = xvec( mask );
            vol_form = vol_form( mask );     
        end
        
        % estimate of L1 by integrating volume form over the domain using
        % the trapezoid rule
        L(1) = diff( xvec ) * ( vol_form(1:end-1) + vol_form(2:end) ) / 2;
        
    case 2
        % get the voxel grid dimensions after resolution increase
        dx = diff( xvals{1} );
        dx = dx(1);
        dy = diff( xvals{2} );
        dy = dy(1);
        
        % short cuts for the metric entries
        g_xx = g(:, :, 1, 1 );
        g_yy = g(:, :, 2, 2 );
        g_xy = g(:, :, 1, 2 );
        
        clear g
        
        %%% calculate LKC 2
        % get the volume form, max introduced for stability, since at the
        % boundaries there might be tiny negative numbers due to numerical
        % calculations
        vol_form = sqrt( max( g_xx .* g_yy - g_xy.^2, 0 ) );
        
        % restrict vol_form to mask
        if mask_opt(2) == 1
            vol_form = vol_form( mask );
        end
        
        % integate volume form over the domain. It assumes that each voxel
        % has the same volume dx*dy and simple midpoint integration is used
        % note that we also use weights to give an appropriate quotient to
        % boundary voxels, if resAdd > 0, which takes into account that the
        % volume of the boundary voxels needs to be halved, or quartered, etc
        L(2) = sum( vol_form(:) .* weights(:) ) * dx * dy;
            
        %%% calculate LKC 1
        % find x shift boundary voxels, i.e. horizontal boundary parts, and
        % integrate using trapozoid rule.
        % Note that we later need to remove half of the end points value,
        % of each line segment which are contained in the x and the y shift
        % boundary. They will be count double otherwise.
        xbdry = bdry_voxels( mask, "x" );      
        L(1)  = sum( sqrt( g_xx( xbdry ) ) ) * dx;
        
        % find y shift boundary voxels, i.e. vertical boundary parts, and
        % integrate using trapozoid rule.
        ybdry = bdry_voxels( mask, "y" );      
        L(1)  = L(1) + sum( sqrt( g_yy( ybdry ) ) ) * dy;
        
        % remove double counted voxels at end of line segments and divide
        % length of boundary by 2, since that is what LKC1 is.
        xybdry = ybdry + xbdry == 2;
        L(1) = ( L(1) - sum( sqrt( g_xx( xybdry ) ) ) * dx / 2 ...
                      - sum( sqrt( g_yy( xybdry ) ) ) * dy / 2 ) / 2;
        
    case 3
        % get the voxel grid dimensions after resolution increase
        dx = diff( xvals{1} );
        dx = dx(1);
        dy = diff( xvals{2} );
        dy = dy(1);
        dz = diff( xvals{3} );
        dz = dz(1);
        
        % short cuts for the metric entries
        g_xx = g(:, :, :, 1, 1 );
        g_yy = g(:, :, :, 2, 2 );
        g_zz = g(:, :, :, 3, 3 );
        g_xy = g(:, :, :, 1, 2 );
        g_xz = g(:, :, :, 1, 3 );
        g_yz = g(:, :, :, 2, 3 );
        
        %%% calculate LKC 3
        % get the volume form, i.e. sqrt(det g), max introduced for
        % stability
        vol_form = sqrt( max(   g_xx.*g_yy.*g_zz...
                              + g_xy.*g_yz.*g_xz...
                              + g_xz.*g_xy.*g_yz...
                              - g_xz.^2.*g_yy...
                              - g_xy.^2.*g_zz...
                              - g_xx.*g_yz.^2, 0 ) );

        % restrict vol_form to mask
        if mask_opt(2) == 1 
            vol_form = vol_form( mask );
        end
        
        % integate volume form over the domain assuming each voxel having
        % the same volume dxdydz. Simple midpoint integration is used.
        L(3) = sum( vol_form(:) * weights(:) ) * dx * dy * dz;                  

        %%% calculate LKC 2
        % find faces having constant z value and integrate using simple
        % midpoint rule.
        bdry = bdry_voxels( mask, "xy" );
        L(2) = sum( sqrt( max( g_xx( bdry ) .* g_yy( bdry )...
                              - g_xy( bdry ).^2, 0 ) ) ) * dx * dy / 2;
                          
        bdry = bdry_voxels( mask, "xz" );
        L(2) = L(2) + sum( sqrt( max( g_xx( bdry ) .* g_zz( bdry )...
                              - g_xz( bdry ).^2, 0 ) ) * dx * dz / 2 );

        bdry = bdry_voxels( mask, "yz" );
        L(2) = L(2) + sum( sqrt( max( g_yy( bdry ) .* g_zz( bdry )...
                              - g_yz( bdry ).^2, 0 ) ) * dy * dz / 2 );
        
        %%% calculate LKC 1
        % work in progress
end
%%%%%% END estimate the LKCs in different dimensions
           
return