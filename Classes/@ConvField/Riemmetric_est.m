function g = riemmetric_est( cfield, dcfield )
% riemmetric_est( cfield, dcfield ) estimates the induced Riemannian metric
% from a convolution field and its derivative.
% This function uses the exact formula from taking the derivatives of the
% kernels and computing the corresponding convolutions, see [1, Theorem ???].
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  cfield   an object of class ConvField of derivtype 0, fiberD = 1 and
%           fibersize > 1.
%  dcfield  an object of class ConvField of derivtype 1, fiberD = 1 and
%           fibersize > 1.
%
%--------------------------------------------------------------------------
% OUTPUT
%   g    an object of class Field representing the Riemannian metric
%        induced by the convolution field.
% -------------------------------------------------------------------------
% DEVELOPER TODOs:
% - write meaningful examples
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow
%--------------------------------------------------------------------------

%% Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%% Check the cfield input
if cfield.fiberD ~= 1
    error( "The cfield input must be multiple observations of a scalar field." )
end

if cfield.fibersize == 1
    error( "The cfield input must have more than one observation of a scalar field." )
end

if cfield.derivtype ~= 0 
    error( "The cfield must be of derivtype 0." )
end


%%% Check the dcfield input
if dcfield.fiberD ~= cfield.fiberD + 1 && dcfield.fiberD ~= 1 
    error( "The dcfield input must be multiple observations of a scalar field." )
end

if dcfield.fibersize(1) == 1
    error( "The dcfield input must have more than one observation of a scalar field." )
end

if dcfield.derivtype ~= 1
    error( "The dcfield must be of derivtype 1." )
end

%%% Check that cfield and dcfield are compatible
if cfield.D ~= dcfield.D
    error( "cfield and dcfield need to be defined on a domain with the same dimension." )
end

if cfield.masksize ~= dcfield.masksize
    error( "cfield and dcfield need to be defined on the same domain." )
end


% Check that method is implemented for dimension D
if cfield.D > 3
    error( 'D must be < 4. Higher dimensional domains have not been implemented')
end

%%% Get constants
% Dimension of the domain
D = cfield.D;

% Get number of subjects/samples
nsubj = cfield.fibersize;

% Get variable domain counter
index  = repmat( {':'}, 1, D );

%% Main function
%--------------------------------------------------------------------------

% Allocate output Field object for the Riemannian metric by initializing it
% with the mask from cfield.
g = Field( cfield.mask );
% Insert xvals values into the Field and initialise the field property with
% the correct dimensions
g.xvals = cfield.xvals;
g.field = zeros( [ cfield.masksize D D ] );

%%%%%% BEGIN compute the Riemannian metric
switch D
    case 1
        % Get the estimates of the variances and covariances required to
        % compute the Riemannian metric
        VY    = var( cfield.field,  0, D+1 );
        VdxY  = var( dcfield.field, 0, D+1 );
        CYdxY = sum( ( dcfield.field - mean( dcfield.field, D+1 ) )...
                                .* cfield.field, D+1 ) / ( nsubj - 1 );
        
        % Get the Riemannian metric at each point
        g.field( index{:}, 1 ) = ( VdxY .* VY - CYdxY.^2 ) ./ VY.^2;
        
    case 2
        % Rename the partial derivatives of the convolution field
        convY  = cfield.field;
        convYx = squeeze( dcfield.field( index{:}, :, 1 ) );
        convYy = squeeze( dcfield.field( index{:}, :, 2 ) );
        clear cfield dcfield
        
        % Get the estimates of the variances and covariances required to
        % compute the Riemannian metric
        VY   = var( convY,  0, D+1 ); 
        VdxY = var( convYx, 0, D+1 );
        VdyY = var( convYy, 0, D+1 );
        CdxYdyY = sum( ( convYy - mean( convYy, D+1 ) ) .* ...
                         convYx, D+1 ) / ( nsubj - 1 );
        CYdxY = sum( ( convYx - mean( convYx, D+1 ) ) .* ...
                       convY, D+1 ) / ( nsubj - 1 );
        CYdyY = sum( ( convYy - mean( convYy, D+1 ) ) .* ...
                       convY, D+1 ) / ( nsubj - 1 );
                 
        % Entries of Riemanian metric
        g_xx = -CYdxY.^2 ./ VY.^2 + VdxY ./ VY;
        g_yy = -CYdyY.^2 ./ VY.^2 + VdyY ./ VY;
        g_xy = -CYdyY .* CYdxY ./ VY.^2 + CdxYdyY ./ VY;
        
        % Prepare output
        g.field( index{:}, 1, 1 ) = g_xx;
        g.field( index{:}, 1, 2 ) = g_xy;
        g.field( index{:}, 2, 1 ) = g_xy;
        g.field( index{:}, 2, 2 ) = g_yy;
        
    case 3
        % Rename the partial derivatives of the convolution field
        convY  = cfield.field;
        convYx = squeeze( dcfield.field( index{:}, :, 1 ) );
        convYy = squeeze( dcfield.field( index{:}, :, 2 ) );
        convYz = squeeze( dcfield.field( index{:}, :, 3 ) );
        clear cfield dcfield
        
        % Get the estimates of the covariances
        VY   = var( convY,  0, D+1 );
        VdxY = var( convYx, 0, D+1 );
        VdyY = var( convYy, 0, D+1 );
        VdzY = var( convYz, 0, D+1 );
        
        CdxYdyY = sum( ( convYy - mean( convYy, D+1 ) ) .* ...
                         convYx, D+1 ) / (nsubj-1);
        CdxYdzY = sum( ( convYx - mean( convYx, D+1 ) ) .* ...
                         convYz, D+1 ) / (nsubj-1);
        CdyYdzY = sum( ( convYy - mean( convYy, D+1 ) ) .* ...
                         convYz, D+1 ) / (nsubj-1);

        CYdxY = sum( ( convYx - mean( convYx, D+1 ) ) .* ...
                       convY, D+1 ) / (nsubj-1);
        CYdyY = sum( ( convYy - mean( convYy, D+1 ) ) .* ...
                       convY, D+1 ) / (nsubj-1);
        CYdzY = sum( ( convYz - mean( convYz, D+1 ) ) .* ...
                       convY, D+1 ) / (nsubj-1);
                 
        % Entries of Riemanian metric
        g_xx = ( -CYdxY.^2 + VdxY .* VY ) ./ VY.^2;
        g_yy = ( -CYdyY.^2 + VdyY .* VY ) ./ VY.^2;
        g_zz = ( -CYdzY.^2 + VdzY .* VY ) ./ VY.^2;
        g_xy = ( -CYdyY .* CYdxY + CdxYdyY .* VY ) ./ VY.^2;
        g_xz = ( -CYdzY .* CYdxY + CdxYdzY .* VY ) ./ VY.^2;
        g_yz = ( -CYdzY .* CYdyY + CdyYdzY .* VY ) ./ VY.^2;
        
        % Prepare output
        g.field( index{:}, 1, 1 ) = g_xx;
        g.field( index{:}, 1, 2 ) = g_xy;
        g.field( index{:}, 1, 3 ) = g_xz;
        g.field( index{:}, 2, 1 ) = g_xy;
        g.field( index{:}, 2, 2 ) = g_yy;
        g.field( index{:}, 2, 3 ) = g_yz;
        g.field( index{:}, 3, 1 ) = g_xz;
        g.field( index{:}, 3, 2 ) = g_yz;
        g.field( index{:}, 3, 3 ) = g_zz;
        
end
%%%%%% END compute the Riemannian metric

% Remove NaNs and replace with zero
g.field = nan2zero( g.field );

return