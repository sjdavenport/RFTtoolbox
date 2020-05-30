function [ Lambda_est, xvals ] = Lambda_conv_est( lat_data, Kernel, resAdd,...
                                                  shift_bdry )
% This function calculates an estimate of the induced Riemannian metric
% from a convolution field or also known as Lambda matrix in neuroimaging
% by using the exact formula from taking the derivatives of the kernels.
%
% It is heavily using the function convfield() and allows for the same
% input to make it compatible.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   lat_data data array T_1 x ... x T_D x N. Last index enumerates the
%            samples. Note that N > 1 is required!
%   Kernel   array 1x1 or 1xD containing the FWHM for different directions
%            for smoothing with a Gaussian kernel, if numeric an isotropic
%            kernel is assumed.
% Optional
%   resAdd   integer denoting the amount of voxels padded between existing
%            voxels to increase resolution
%   shift_bdry integer denoting the amount of voxels padded between existing
%              voxels to increase resolution
%--------------------------------------------------------------------------
% OUTPUT
%   Lambda_est  an array of size hrT_1 x ... x hrT_D x D x D. Containing
%               the D x D Lambda matrix/Riemannian metric at each voxel
%               from the chosen increased resolution defined by resAdd.
%   xvals       the grid values of the high resolution grid.
% -------------------------------------------------------------------------
% DEVELOPER TODOs:
%   - add same input as convfield as soon as convfield is finalized
%--------------------------------------------------------------------------
% AUTHORS: Fabian Telschow
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
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Check input and get important constants from the mandatory input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get constants from the mandatory input
% size of the domain
s_lat_data = size( lat_data );

% dimension of the domain
D = length( s_lat_data( 1:end-1 ) );

% get number of subjects/samples
nsubj = s_lat_data( D + 1 );

% get variable domain counter
index  = repmat( {':'}, 1, D );

% check that method is implemented for dimension D
if D > 3
    error( 'D must be < 4. Higher dimensional domains have not been implemented')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% add/check optional values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist( 'resAdd', 'var' )
   % default number of resolution increasing voxels between observed voxels
   resAdd = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% main function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% get the convolution fields and their derivatives
[ convY, xvals ] = convfield( lat_data, Kernel, D, resAdd, shift_bdry );
DconvY           = convfield( lat_data, Kernel, D, resAdd, 1, shift_bdry );

sY = size( convY );

%%%% allocate output the entries of the Riemannian metric
Lambda_est = NaN * ones( [ sY(1:end-1)  D D ] );

%%%%%%%% BEGIN compute the Riemannian metric
switch D
    case 1
        % Get the estimates of the variances and covariances required to
        % compute the Riemannian metric/Lambda matrix
        VY    = var( convY,  0, D+1 );
        VdxY  = var( DconvY, 0, D+1 );
        CYdxY = sum( ( DconvY - mean( DconvY, D+1 ) ) .* convY, D+1 ) ...
                                                            / ( nsubj - 1 );
        
        % get the Riemannian metric at each point
        Lambda_est( index{:}, 1 ) = ( VdxY .* VY - CYdxY.^2 ) ./ VY.^2;
    case 2
        % rename the partial derivatives of the convolution field
        convYx = squeeze( DconvY( 1, index{:}, : ) );
        convYy = squeeze( DconvY( 2, index{:}, : ) );
        clear DconvY
        
        % Get the estimates of the variances and covariances required to
        % compute the Riemannian metric/Lambda matrix
        VY   = var( convY,  0, D+1 ); 
        VdxY = var( convYx, 0, D+1 );
        VdyY = var( convYy, 0, D+1 );
        CdxYdyY = sum( ( convYy - mean( convYy, D+1 ) ) .* ...
                         convYx, D+1 ) / ( nsubj - 1 );
        CYdxY = sum( ( convYx - mean( convYx, D+1 ) ) .* ...
                       convY, D+1 ) / ( nsubj - 1 );
        CYdyY = sum( ( convYy - mean( convYy, D+1 ) ) .* ...
                       convY, D+1 ) / ( nsubj - 1 );
                 
        % entries of Riemanian metric
        g_xx = -CYdxY.^2 ./ VY.^2 + VdxY ./ VY;
        g_yy = -CYdyY.^2 ./ VY.^2 + VdyY ./ VY;
        g_xy = -CYdyY .* CYdxY ./ VY.^2 + CdxYdyY ./ VY;
        
        % prepare output
        Lambda_est( index{:}, 1, 1 ) = g_xx;
        Lambda_est( index{:}, 1, 2 ) = g_xy;
        Lambda_est( index{:}, 2, 1 ) = g_xy;
        Lambda_est( index{:}, 2, 2 ) = g_yy;
    case 3
        % rename the partial derivatives of the convolution field
        convYx = squeeze( DconvY( 1, index{:}, : ) );
        convYy = squeeze( DconvY( 2, index{:}, : ) );
        convYz = squeeze( DconvY( 3, index{:}, : ) );
        clear DconvY
        
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
        g_xx = ( -CYdxY.^2 + VdxY .* VY ) ./ VY.^2;
        g_yy = ( -CYdyY.^2 + VdyY .* VY ) ./ VY.^2;
        g_zz = ( -CYdzY.^2 + VdzY .* VY ) ./ VY.^2;
        g_xy = ( -CYdyY .* CYdxY + CdxYdyY .* VY ) ./ VY.^2;
        g_xz = ( -CYdzY .* CYdxY + CdxYdzY .* VY ) ./ VY.^2;
        g_yz = ( -CYdzY .* CYdyY + CdyYdzY .* VY ) ./ VY.^2;
        
        % prepare output
        Lambda_est( index{:}, 1, 1 ) = g_xx;
        Lambda_est( index{:}, 1, 2 ) = g_xy;
        Lambda_est( index{:}, 1, 3 ) = g_xz;
        Lambda_est( index{:}, 2, 1 ) = g_xy;
        Lambda_est( index{:}, 2, 2 ) = g_yy;
        Lambda_est( index{:}, 2, 3 ) = g_yz;
        Lambda_est( index{:}, 3, 1 ) = g_xz;
        Lambda_est( index{:}, 3, 2 ) = g_yz;
        Lambda_est( index{:}, 3, 3 ) = g_zz;
end
%%%%%%%% END compute the Riemannian metric
return