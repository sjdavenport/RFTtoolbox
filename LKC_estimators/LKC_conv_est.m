function LKC = LKC_conv_est( lat_data, mask, Kernel, resAdd, mask_opt,...
                             enlarge, version )
% LKC_CONV_EST( lat_data, mask, Kernel, resAdd, mask_opt, enlarge, version )
% estimates the Lipschitz Killing curvatures for a convolution field
% derived from a general kernel. Corresponding samples of such fields can
% be simulated using convfield.m which has the same input variables.
% If 'version' is choosing to be 'analytical' the analytic derivatives in
% form of convolutions with the derivative of the kernel are used. If
% version is 'numerical', numerical approximation of the derivative are
% used, which are precise approximations of the true derivative, since
% convolution fields can be compute for any location.
%
% Currently, only 1D and 2D provide estimates for all LKCs.
% 3D only allows for estimation of L2 and L3.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   lat_data  data array T_1 x ... x T_D x N. Last index enumerates the
%             samples. Note that N > 1 is required!
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
%   version    string indicating which estimator for the Lambda
%              matrix/Riemannian metric is used. Options are "analytical"
%              and "numerical". Default "analytical".
%              Note that "numerical" might be faster.
%--------------------------------------------------------------------------
% OUTPUT
%   LKC  structure containing fields:
%        - hatL: D x 1 vector of estimates of LKC for the sample Y. It
%                is the average of hatL1.
%        - L0:   integer containing the Euler characteristic of the mask
%                equivalently the zeroth LKC of the random fields.
%        - geom: structure containing geometric quantities of the
%                induced metric of the random field.
%--------------------------------------------------------------------------
% DEVELOPER TODOs:
%   - add full 3D estimation
% -------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow
%--------------------------------------------------------------------------


%% %-----------------------------------------------------------------------
%  check mandatory input and get important constants
%--------------------------------------------------------------------------
% size of the domain
s_lat_data = size( lat_data );

% get size of the mask
sM = size( mask );

% design question, do we want to force the user to use logicals?
mask = logical( mask );

% dimension of the domain, since matlab is not consistent for D<1, we need
% to catch this case 
if sM( 2 ) == 1 && length( sM ) == 2
    D = 1;
else
    D = length( sM );
end

% check whether more than one sample is provided, else reject the input.
if length( s_lat_data ) <= D || ( D == 1 && s_lat_data(2) == 1 )
    error( "At least 2 samples of a random field are needed!" );
elseif length( s_lat_data ) > D + 1
    error( "Y needs to have only one dimension more than mask!" );
else
    % get number of subjects/samples
    nsubj = s_lat_data( D + 1 );
end

% check validity of mask input
if ~all( sM == s_lat_data( 1:end-1 ) ) && ~all( sM == s_lat_data ) && ...
   ~( D == 1 && sM(1) == s_lat_data(1) )
   error( 'The mask needs to have the same size as Y.\n%s',...
          'Note that realisations are stored as columns.' )
end

% check that method is implemented for dimension D
if D > 3
    error( 'D must be < 4. Higher dimensional domains have not been implemented')
end

%% %-----------------------------------------------------------------------
%  add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'resAdd', 'var' )
   % default number of resolution increasing voxels between observed voxels
   resAdd = 1;
end

if ~exist( 'version', 'var' )
   % default method for Lambda matrix estimation
   version = "analytical";
end

if ~exist( 'enlarge', 'var' )
   % default method for Lambda matrix estimation
   enlarge = ceil( resAdd / 2 );
end

if ~exist( 'mask_opt', 'var' )
   % default method for mask_opt, which controls when the mask is applied
   % prior/after application of smoothing using convfield.m
   % Currently we recommend using it prior and after, yet this might depend
   % on your application!
   mask_opt = [ 1 1 ];
end

%% %-----------------------------------------------------------------------
%  main function
%--------------------------------------------------------------------------
%%% allocate variables
% allocate vector for Lipschitz Killing curvature
L = NaN * ones( [ 1 D ] );
% structure to output the different computed fields for debugging purposes
geom = struct();

%%% manipulate the mask and mask the data
% mask the lattice data, if opted for
if mask_opt(1) == 1
    lat_data = repmat( mask, [ ones( [ 1 D ] ), nsubj ] ) .* lat_data;
end

% get mask on higher resolution and the weights of each voxel for the
% volume computation.
if resAdd ~= 0
    [ mask, weights ] = mask_highres( mask, resAdd, enlarge );
    % reduce weights matrix only to active voxels for speed up
    weights = weights( mask );
else
    weights = ones( [ sum( mask(:) ), 1 ] );
end

%%% get the Riemannian metric/Lambda matrix from the data
if strcmp( version, "analytical")
    [ g, xvals ] = Lambda_conv_est( lat_data, Kernel, resAdd, enlarge );
elseif strcmp( version, "numerical")
    [ g, xvals ] = Lambda_numeric_est( lat_data, Kernel, resAdd, enlarge );
else
    % output error message, if not a valid method is chosen
    error( strcat( "Choose a valid method for Lambda/Riemannian metric", ...
           " estimation.\n%s" ),...
           "Options are 'analytical' and 'numerical'." )
end


%%%%%% BEGIN estimate the LKCs in different dimensions
%%% compute 0th LKC
L0 = EulerChar( mask, 0.5, D );

%%% compute LKCs for 0 < d < D
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
        
        %%% fill the output structure
        geom.vol_form    = vol_form;
        geom.riem_metric = g;
        
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
        
        % save g to the output structure and clear
        geom.riem_metric = g;
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
                        
        %%% fill the output structure
        geom.vol_form = vol_form;
        
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
        
        % save g to the output structure and clear
        geom.riem_metric = g;
        clear g
        
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


%% %-----------------------------------------------------------------------
%  prepare output structure
%--------------------------------------------------------------------------
% summarize output
LKC  = struct( 'hatL', L, 'L0', L0, 'geomQuants',...
               geom );
           
return