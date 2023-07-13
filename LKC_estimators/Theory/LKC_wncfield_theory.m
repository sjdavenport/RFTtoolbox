function [LKC, g] = LKC_wncfield_theory( mask, params, mask2 )
% LKC_WNCFIELD_THEORY( mask, params ) computes
% theoretical Lipschitz Killing curvatures for a convolution field
% derived from a seperable kernel with underlying discrete independent
% Gaussian white noise process on a lattice.
% It uses the fact that the voxels are independent and the variance and
% mean are known.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  mask     a logical array of dimension T_1 x...x T_D.
%  params   an object of class ConvFieldParams
%  mask2    optional second mask to allow for more non-stationarity
%
%--------------------------------------------------------------------------
% OUTPUT
%   L  1 x D array of theoretical LKCs computed assuming independence
%      of voxels mean zero and variance 1
%--------------------------------------------------------------------------
% EXAMPLES
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow
%--------------------------------------------------------------------------

%% check mandatory input and get important constants
%--------------------------------------------------------------------------

% Get size of the mask
sM = size( mask );
% make mask logical
mask = logical( mask );

if ~exist( 'mask2', 'var' )
   % default number of bootstrap replicates
   mask2 = mask;
end

% Dimension of the domain, since matlab is not consistent for D<1, we need
% to catch this case 
if sM( 2 ) == 1 && length( sM ) == 2
    D = 1;
else
    D = length( sM );
end

% Check that method is implemented for dimension D
if D > 3
    error( strcat( 'D must be < 4. Higher dimensional domains have',...
                   'not been implemented' ) );
end

% Save the parameter from params in a small file
Kernel     = params.kernel;
resadd     = params.resadd;
enlarge    = params.enlarge;
lat_masked = params.lat_masked;


%% Add/check optional input
%--------------------------------------------------------------------------

%% Main function
%--------------------------------------------------------------------------

% Get dx from resadd
dx = 1 / ( resadd + 1 );

% Get mask on higher resolution and the weights of each voxel for the
% volume computation.
if resadd ~= 0
    [ mask, ~ ] = mask_highres( mask, resadd, enlarge );
    [ mask2, ~ ] = mask_highres( mask2, resadd, enlarge );
end

% Get the size of the resolution increased domain
Dimhr = size( mask );


%%%%%% BEGIN get the Riemannian metric/Lambda matrix assuming that the
%%%%%% voxels are independent Gaussian white noise
% Preallocate the Riemannian metric
g = zeros( [ Dimhr D D ] );

% Create index to fill the original mask at the correct voxels of the high
% resolution mask
index = cell( [ 1 D ] );
for d = 1:D
    index{d} = ( enlarge + 1 ):( resadd + 1 ):( Dimhr(d) - enlarge );
end

% Fill the resolution increased mask with the values from the original mask
% at the correct locations
onesField = zeros( Dimhr );
onesField( index{:} ) = 1;

% Mask the field, if opted for
if lat_masked
    onesField = mask .* onesField;
end

%%% Get the entries of the Riemannian metric
switch D
    case 1
        % Get the theoretical variance of the field and the variance of
        % derivatives
        VY    = fconv( onesField, @(x) Kernel.kernel{1}(x).^2, D,...
                        Kernel.truncation, dx );
        VdxY  = fconv( onesField, @(x) Kernel.dkernel{1}(x).^2, D,...
                        Kernel.truncation, dx );
        CYdxY = fconv( onesField, @(x) Kernel.kernel{1}(x) .* ...
                                        Kernel.dkernel{1}(x), D,...
                        min( Kernel.truncation,  Kernel.dtruncation ), dx );
        
        % Get the volume form
        g = ( VdxY .* VY - CYdxY.^2 ) ./ VY.^2;
        
    case 2
        % Get the partial derivatives of the Kernel
        grad_Kernel = Gradient( Kernel );
        dxKernel = grad_Kernel{1};
        dyKernel = grad_Kernel{2};
        clear grad_Kernel
        
        % Get the theoretical variance of the field and the variance of
        % derivatives
        VY      = fconv( onesField,...
                          { @(x) Kernel.kernel{1}(x).^2,...
                            @(x) Kernel.kernel{2}(x).^2 }, D,...
                          Kernel.truncation, dx );
        VdxY    = fconv( onesField,...
                          { @(x) dxKernel.kernel{1}(x).^2,...
                            @(x) dxKernel.kernel{2}(x).^2 }, D,...
                          Kernel.truncation, dx );
        VdyY    = fconv( onesField,...
                          { @(x) dyKernel.kernel{1}(x).^2,...
                            @(x) dyKernel.kernel{2}(x).^2 }, D,...
                          Kernel.truncation, dx );
        CYdyY   = fconv( onesField,...
                          { @(x) Kernel.kernel{1}(x) .* dyKernel.kernel{1}(x),...
                            @(x) Kernel.kernel{2}(x) .* dyKernel.kernel{2}(x) },...
                            D, Kernel.truncation, dx );
        CYdxY   = fconv( onesField,...
                          { @(x) Kernel.kernel{1}(x) .* dxKernel.kernel{1}(x),...
                            @(x) Kernel.kernel{2}(x) .* dxKernel.kernel{2}(x) },...
                            D, Kernel.truncation, dx );
        CdxYdyY = fconv( onesField,...
                          { @(x) dxKernel.kernel{1}(x) .* dyKernel.dkernel{1}(x),...
                            @(x) dxKernel.kernel{2}(x) .* dyKernel.dkernel{2}(x) },...
                            D, Kernel.truncation, dx );

        % Entries of riemanian metric
        g( :, :, 1, 1 ) = -CYdxY.^2 ./ VY.^2 + VdxY ./ VY;
        g( :, :, 2, 2 ) = -CYdyY.^2 ./ VY.^2 + VdyY ./ VY;
        g( :, :, 1, 2 ) = -CYdyY .* CYdxY ./ VY.^2 + CdxYdyY ./ VY;
        g( :, :, 1, 2 ) = g( :, :, 2, 1 );

     case 3
        % Get the partial derivatives of the Kernel
        grad_Kernel = Gradient( Kernel );
        dxKernel = grad_Kernel{1};
        dyKernel = grad_Kernel{2};
        dzKernel = grad_Kernel{3};
        clear grad_Kernel
        
        % Get the estimates of the covariances        
        VY   = fconv( onesField,...
                       { @(x) Kernel.kernel{1}(x).^2,...
                         @(x) Kernel.kernel{2}(x).^2,...
                         @(x) Kernel.kernel{3}(x).^2 }, D,...
                       Kernel.truncation, dx );
        VdxY = fconv( onesField,...
                       { @(x) dxKernel.kernel{1}(x).^2,...
                         @(x) dxKernel.kernel{2}(x).^2,...
                         @(x) dxKernel.kernel{3}(x).^2 }, D,...
                       Kernel.truncation, dx );
        VdyY = fconv( onesField,...
                       { @(x) dyKernel.kernel{1}(x).^2,...
                         @(x) dyKernel.kernel{2}(x).^2,...
                         @(x) dyKernel.kernel{3}(x).^2 }, D,...
                       Kernel.truncation, dx );
        VdzY = fconv( onesField,...
                       { @(x) dzKernel.kernel{1}(x).^2,...
                         @(x) dzKernel.kernel{2}(x).^2,...
                         @(x) dzKernel.kernel{3}(x).^2 }, D,...
                       Kernel.truncation, dx );
        
        CdxYdyY = fconv( onesField,...
                          { @(x) dxKernel.kernel{1}(x) .* dyKernel.kernel{1}(x),...
                            @(x) dxKernel.kernel{2}(x) .* dyKernel.kernel{2}(x),...
                            @(x) dxKernel.kernel{3}(x) .* dyKernel.kernel{3}(x) },...
                            D, Kernel.truncation, dx );

        CdxYdzY = fconv( onesField,...
                          { @(x) dxKernel.kernel{1}(x) .* dzKernel.kernel{1}(x),...
                            @(x) dxKernel.kernel{2}(x) .* dzKernel.kernel{2}(x),...
                            @(x) dxKernel.kernel{3}(x) .* dzKernel.kernel{3}(x) },...
                            D, Kernel.truncation, dx );
        CdyYdzY = fconv( onesField,...
                          { @(x) dzKernel.kernel{1}(x) .* dyKernel.kernel{1}(x),...
                            @(x) dzKernel.kernel{2}(x) .* dyKernel.kernel{2}(x),...
                            @(x) dzKernel.kernel{3}(x) .* dyKernel.kernel{3}(x) },...
                            D, Kernel.truncation, dx );
                        
        CYdxY = fconv( onesField,...
                        { @(x) dxKernel.kernel{1}(x) .* Kernel.kernel{1}(x),...
                          @(x) dxKernel.kernel{2}(x) .* Kernel.kernel{2}(x),...
                          @(x) dxKernel.kernel{3}(x) .* Kernel.kernel{3}(x) },...
                          D, Kernel.truncation, dx );
        CYdyY = fconv( onesField,...
                        { @(x) dyKernel.kernel{1}(x) .* Kernel.kernel{1}(x),...
                          @(x) dyKernel.kernel{2}(x) .* Kernel.kernel{2}(x),...
                          @(x) dyKernel.kernel{3}(x) .* Kernel.kernel{3}(x) },...
                          D, Kernel.truncation, dx );
        CYdzY = fconv( onesField,...
                        { @(x) dzKernel.kernel{1}(x) .* Kernel.kernel{1}(x),...
                          @(x) dzKernel.kernel{2}(x) .* Kernel.kernel{2}(x),...
                          @(x) dzKernel.kernel{3}(x) .* Kernel.kernel{3}(x) },...
                          D, Kernel.truncation, dx );
                 
        % Entries of riemanian metric/ Lambda matrix from neuroimaging
        g( :, :, :, 1, 1 ) = ( -CYdxY.^2 + VdxY .* VY ) ./ VY.^2;
        g( :, :, :, 2, 2 ) = ( -CYdyY.^2 + VdyY .* VY ) ./ VY.^2;
        g( :, :, :, 3, 3 ) = ( -CYdzY.^2 + VdzY .* VY ) ./ VY.^2;
        g( :, :, :, 1, 2 ) = ( -CYdyY .* CYdxY + CdxYdyY .* VY ) ./ VY.^2;
        g( :, :, :, 2, 1 ) = g( :, :, :, 1, 2 );
        g( :, :, :, 1, 3 ) = ( -CYdzY .* CYdxY + CdxYdzY .* VY ) ./ VY.^2;
        g( :, :, :, 3, 1 ) = g( :, :, :, 1, 3 );
        g( :, :, :, 2, 3 ) = ( -CYdzY .* CYdyY + CdyYdzY .* VY ) ./ VY.^2;
        g( :, :, :, 3, 2 ) = g( :, :, :, 2, 3 );
        
end

% Change nan to zero because there is a division by zero if masking
g = nan2zero( g );

% Setting up the default xvals_vecs
xvals  = cell( 1, D );

for d = 1:D
    xvals{d} = ( ( 1 - enlarge*dx ):dx:( sM(d) + enlarge*dx ) )...
                + Kernel.adjust(d);
end

% Create a voxmanifold with the computed metric
tmp = g;
g   = Field( mask2 );
g.field = tmp;
g.xvals = xvals;
g = Mask(g);
voxmfd  = VoxManifold( g );
           
LKC = LKC_est( voxmfd );

return