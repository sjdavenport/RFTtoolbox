function LKC = LKC_conv_est( lat_data, mask, Kernel, resadd, mask_lat,...
                             enlarge, version )
% LKC_CONV_EST( lat_data, mask, Kernel, resadd, mask_lat, enlarge, version )
% estimates the Lipschitz Killing curvatures for a convolution field
% derived from a general kernel. Corresponding samples of such fields can
% be simulated using convfield.m which has the same input variables.
% If 'version' is choosen to be 'analytical' the analytic derivatives in
% form of convolutions with the derivatives of the kernel are used. If
% version is 'numerical', numerical approximation of the derivatives are
% used, which are precise approximations of the true derivative, since
% convolution fields can be computed for any location.
%
% Theoretical values for the LKCs of such processes under the assumption
% that lat_data has independent mean zero unit variance voxels can be
% obtained from the function LKC_wncfield_theory.m
%
% While 1D and 2D allow for general non-isotropic fields, the estimate in
% 3D assumes 'local isotropy' for L1. This might be fixed in the future.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  lat_data  data array T_1 x ... x T_D x N. Last index enumerates the
%            samples. Note that N > 1 is required!
%  mask      a logical array of dimension T_1 x...x T_D indicating which
%            voxels belong to the domain of the random field.
%  Kernel    either an object of class SepKernel or a numeric.
%            If class SepKernel all fields must be specified.
%
%            If Kernel is numeric, the convolution field is generated by 
%            smoothing with an isotropic Gaussian kernel with FWHM = Kernel.
%            Truncation and adjust_kernel are set to be default values.
% Optional
%  resadd     integer denoting the amount of voxels padded between 
%             existing voxels to increase resolution. Default 1.
%  mask_lat   an logical. If 1 the mask is applied to the lat_data prior
%             to any other calculation. Default 0.
%             Note that if 0 the mask is only applied after convolutions
%             are performed.              
%  enlarge    an integer denoting the amount of voxels the resolution
%             increased mask is enlarged by dilation. Note that for
%             unbiased estimation of LKCs resadd needs to be an odd number
%             and enlarge needs to be set to the default value.
%             Default ceil( resadd / 2 )
%  version    string indicating which estimator for the Lambda
%             matrix/Riemannian metric is used. Options are "analytical"
%             and "numerical". Default "analytical".
%             Note that "numerical" might be faster.
%--------------------------------------------------------------------------
% OUTPUT
%   LKC  structure containing fields:
%        - hatL: 1 x D vector of estimates of LKC for the sample Y.
%        - L0:   integer containing the Euler characteristic of the mask or
%                equivalently the zeroth LKC of the random fields.
%        - geom: structure containing geometric quantities of the
%                induced metric of the random field.
%--------------------------------------------------------------------------
% DEVELOPER TODOs:
%   - fix local isotropy assumption in L1 by computing the second integral
%--------------------------------------------------------------------------
% EXAMPLES
% %% %% D = 1 
% %% % Field parameters
% T      = 100;
% nsubj  = 120;
% FWHM   = 20;
% pad    = ceil( 4*FWHM2sigma( FWHM ) );
% % Method = "numerical";
% method = "analytical";
% 
% %% Example with recangular mask
% % Generate mask
% mask = pad_vals( true( [ T, 1 ] ), pad, 0 );
% 
% % LKC from continuous theory
% theoryL = LKC_isogauss_theory( FWHM, T  );
% 
% % Generate test data
% lat_data = randn( [ T+2*pad nsubj ] );
% 
% % Closest approximation of the continuous field uses thresholding
% % mask_lat = 0, since otherwise there are boundary effects
% mask_lat = 0;
% LKC1 = LKC_conv_est( lat_data, mask, FWHM, 1, mask_lat );
% LKC3 = LKC_conv_est( lat_data, mask, FWHM, 3, mask_lat );
% LKC5 = LKC_conv_est( lat_data, mask, FWHM, 5, mask_lat );
% 
% % Values are stable accross different resadd increases. Note that resadd
% % should be odd.
% [ theoryL; LKC1.hatL; LKC3.hatL; LKC5.hatL ]'
% 
% % Masking the data shows a boundary effect, validating departure
% % from the theoretical value
% LKC_masked = LKC_conv_est( lat_data, mask, FWHM, 1, 1 );
% [ theoryL; LKC_masked.hatL ]'
% 
% %% %% D = 2 
% %% % Parameters for the field
% T      = 49;
% nsubj  = 100;
% FWHM   = sigma2FWHM(5);
% pad = ceil( 4*FWHM2sigma( FWHM ) );
% 
% %% Example with recangular mask
% % Get mask
% mask = pad_vals( true( [ T T ] ), pad );
% 
% % Get LKC for the theoretical field
% theoryL = LKC_isogauss_theory( FWHM, [ T T ] );
% 
% % Generate test data
% lat_data = randn( [ T+2*pad T+2*pad nsubj ] );
% 
% % Closest approximation of the continuous field uses thresholding
% % mask_lat = 0, since otherwise there are boundary effects
% mask_lat = 0;
% LKC1 = LKC_conv_est( lat_data, mask, FWHM, 1, mask_lat );
% LKC3 = LKC_conv_est( lat_data, mask, FWHM, 3, mask_lat );
% LKC5 = LKC_conv_est( lat_data, mask, FWHM, 5, mask_lat );
% 
% % Values are stable accross different resadd increases. Note that resadd
% % should be odd.
% [ theoryL; LKC1.hatL; LKC3.hatL; LKC5.hatL ]'
% 
% %% Example with complicated mask
% Sig = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
% mask = logical( Sig > 0.02 & Sig < 1.1 );
% figure(1), clf,
% imagesc( mask ), colorbar,
% title("mask")
% clear Sig
% 
% % Get LKC for the theoretical field
% theoryL = LKC_wncfield_theory( mask, FWHM, 3, mask_lat );
% 
% % Generate test data
% lat_data = randn( [ size(mask) nsubj ] );
% 
% % Closest approximation of the continuous field uses thresholding
% % mask_lat = 0, since otherwise there are boundary effects
% mask_lat = 0;
% LKC1 = LKC_conv_est( lat_data, mask, FWHM, 1, mask_lat );
% LKC3 = LKC_conv_est( lat_data, mask, FWHM, 3, mask_lat );
% LKC5 = LKC_conv_est( lat_data, mask, FWHM, 5, mask_lat );
% 
% % Values are stable accross different resadd increases. Note that resadd
% % should be odd.
% [ theoryL.L; LKC1.hatL; LKC3.hatL; LKC5.hatL ]'
% 
% % Masking the data gives different LKCs
% theoryL_masked = LKC_wncfield_theory( mask, FWHM, 3, 1 );
% LKC_masked = LKC_conv_est( lat_data, mask, FWHM, 1, 1 );
% [ theoryL.L; theoryL_masked.L; LKC_masked.hatL ]'
% 
% %% %% D = 3 
% % Parameters for the field
% T      = 20;
% nsubj  = 10;
% FWHM   = sigma2FWHM(1.5);
% pad    = ceil( 4*FWHM2sigma( FWHM ) );
% 
% %% Rectangular domain example
% % Generate rectangular mask with a padded zero collar 
% mask = pad_vals( ones( [ T T T] ), pad );
% 
% % Get theoretical LKC
% theoryL = LKC_wncfield_theory( mask, FWHM, 3, 0 );
% contL = LKC_isogauss_theory( FWHM, [ T T T] );
% 
% % Generate test data
% lat_data = randn( [ T+2*pad T+2*pad T+2*pad nsubj ] );
% 
% % Closest approximation of the continuous field uses thresholding
% % mask_lat = 0, since otherwise there are boundary effects
% mask_lat = 0;
% LKC1 = LKC_conv_est( lat_data, mask, FWHM, 1, mask_lat );
% LKC3 = LKC_conv_est( lat_data, mask, FWHM, 3, mask_lat );
% LKC5 = LKC_conv_est( lat_data, mask, FWHM, 5, mask_lat );
% 
% % Values are stable accross different resadd increases. Note that resadd
% % should be odd.
% [ theoryL.L; contL; LKC1.hatL; LKC3.hatL; LKC5.hatL ]'
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow
%--------------------------------------------------------------------------


%% Check mandatory input and get important constants
%--------------------------------------------------------------------------

% Size of the domain
s_lat_data = size( lat_data );

% Get size of the mask
sM = size( mask );

% Design question, do we want to force the user to use logicals?
mask = logical( mask );

% Dimension of the domain, since matlab is not consistent for D<1, we need
% to catch this case 
if sM( 2 ) == 1 && length( sM ) == 2
    D = 1;
else
    D = length( sM );
end

% Check whether more than one sample is provided, else reject the input.
if length( s_lat_data ) <= D || ( D == 1 && s_lat_data(2) == 1 )
    error( "At least 2 samples of a random field are needed!" );
elseif length( s_lat_data ) > D + 1
    error( "Y needs to have only one dimension more than mask!" );
else
    % Get number of subjects/samples
    nsubj = s_lat_data( D + 1 );
end

% Check validity of mask input
if ~all( sM == s_lat_data( 1:end-1 ) ) && ~all( sM == s_lat_data ) && ...
   ~( D == 1 && sM(1) == s_lat_data(1) )
   error( 'The mask needs to have the same size as Y.\n%s',...
          'Note that realisations are stored as columns.' )
end

% Check that method is implemented for dimension D
if D > 3
    error( strcat( 'D must be < 4. Higher dimensional domains ',...
           'have not been implemented' ) )
end

% If Kernel is numeric use an isotropic Gaussian Kernel else use the
% provided Kernel object
if isnumeric( Kernel ) 
    % Change numerical Kernel input to an object of class SepKernel
    Kernel = SepKernel( D, Kernel );

elseif ~isa( Kernel, 'SepKernel' )
    error( strcat( "The 'Kernel' must be either a numeric or an ",...
                   "object of class SepKernel!" ) );
end

%% add/check optional values
%--------------------------------------------------------------------------

if ~exist( 'resadd', 'var' )
   % Default number of resolution increasing voxels between observed voxels
   resadd = 1;
end

if ~exist( 'enlarge', 'var' )
   % Default method for Lambda matrix estimation
   enlarge = ceil( resadd / 2 );
end

if ~exist( 'mask_lat', 'var' )
   % Default method for mask_lat controlling whether the input lat_data is
   % masked prior to convolving with the kernel.
   mask_lat = 1;
end

if ~exist( 'version', 'var' )
   % Default method for Lambda matrix estimation
   version = "analytical";
end

%% Main function
%--------------------------------------------------------------------------

%%% Allocate variables
% Allocate vector for Lipschitz Killing curvature
L = NaN * ones( [ 1 D ] );
% Structure to output the different computed fields for debugging purposes
geom = struct();

%%% Manipulate the mask and mask the data
% Mask the lattice data, if opted for
if mask_lat
    lat_data = repmat( mask, [ ones( [ 1 D ] ), nsubj ] ) .* lat_data;
end

% Get mask on higher resolution and the weights of each voxel for the
% volume computation.
if resadd ~= 0
    [ mask, weights ] = mask_highres( mask, resadd, enlarge );
    % Reduce weights matrix only to active voxels for speed up
    weights = weights( mask );
else
    weights = ones( [ sum( mask(:) ), 1 ] );
end

%%% Get the Riemannian metric/Lambda matrix from the data
if strcmp( version, "analytical")
    [ g, xvals ] = Lambda_conv_est( lat_data, Kernel, resadd, enlarge );
elseif strcmp( version, "numerical")
    [ g, xvals ] = Lambda_numeric_est( lat_data, Kernel, resadd, enlarge );
else
    % Output error message, if not a valid method is chosen
    error( strcat( "Choose a valid method for Lambda/Riemannian metric", ...
           " estimation.\n%s" ),...
           "Options are 'analytical' and 'numerical'." )
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
        
        %%% Fill the output structure
        geom.vol_form    = vol_form;
        geom.riem_metric = g;
        
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
        
        % Save g to the output structure and clear
        geom.riem_metric = g;
        clear g
        
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
                               
        %%% Fill the output structure
        geom.vol_form = vol_form;
        
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
        
        % Save g to the output structure and clear
        geom.riem_metric = g;
        clear g
        
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
        % Find faces and edges and their integration weights.
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


%% Prepare output structure
%--------------------------------------------------------------------------
% Summarize output
LKC  = struct( 'hatL', L, 'L0', L0, 'geomQuants',...
               geom );
           
return