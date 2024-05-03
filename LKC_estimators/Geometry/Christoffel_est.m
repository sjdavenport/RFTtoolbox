function Gamma = Christoffel_est( field, dfield, d2field )
% Christoffel_est( field, dfield, d2field ) estimates the Christoffel
% symbols of the induced Riemannian metric induced by random twice
% differentiable fields. It uses the analytic formula derived in 
% [1, Theorem ???].
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  field    object of class Field representing observations of the random
%           field with fiberD = 1 and fibersize > 1.
%  dfield   object of class Field representing observations of first
%           derivatives of 'field'.
%  d2field  object of class Field representing observations of second
%           derivatives of 'field'.
%--------------------------------------------------------------------------
% OUTPUT
%  Gamma    object of class Field containing the Christoffel symbols
% -------------------------------------------------------------------------
% DEVELOPER TODOs:
%--------------------------------------------------------------------------
% EXAMPLES
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow
%--------------------------------------------------------------------------

%% Check mandatory input and get important constants
%--------------------------------------------------------------------------

% Check whether the fields have compatible dimensions
if ~iscompatible( field, dfield )
    error( "field and dfield need to be compatible." )
end

% Dimension of the domain
D = field.D;

%% Main function
%--------------------------------------------------------------------------

% Get size of the resolution increased domain
sY = field.masksize;

% Allocate output the entries of the Riemannian metric
Gamma = Field( zeros( [ sY D D D ] ), field.mask, field.xvals );

% Reduce Field objects to arrays
Y       = field.field;
dfield  = dfield.field;
d2field = d2field.field;

% Get variable domain counter
index  = repmat( {':'}, 1, D );

%%%%%% BEGIN compute Christoffel symbols
switch D
    case 1     
        % Rename the partial derivatives of the convolution field
        Yx  = squeeze( dfield( index{:}, :, 1 ) );
        Yxx = squeeze( d2field( index{:}, :, 1, 1 ) );
        
        clear dfield d2field
        
        % Get the estimates of the covariances
        VY   = var( Y,  0, D+1 );
        
        CdYY   = vectcov( Y, Yx, D+1, 1 );
        Cd2YY  = vectcov( Y, Yxx, D+1, 1 );
        CdYdY  = var( Yx,  0, D+1 );
        Cd2YdY = vectcov( Yxx, Yx, D+1, 1 );
        
        Gamma.field( index{:}, 1, 1, 1 ) = ...
                  Cd2YdY ./ VY ...
                - CdYY  .* CdYdY ./ VY.^2 ...
                - Cd2YY .* CdYY  ./ VY.^2 ...
                - CdYY  .* CdYdY ./ VY.^2 ...
                + 2 * CdYY .* CdYY .* CdYY ./  VY.^3;
    case 2
        % Rename the partial derivatives of the convolution field
        Yx = squeeze( dfield( index{:}, :, 1 ) );
        Yy = squeeze( dfield( index{:}, :, 2 ) );
        
        Yxx = squeeze( d2field( index{:}, :, 1, 1 ) );
        Yyy = squeeze( d2field( index{:}, :, 2, 2 ) );
        Yxy = squeeze( d2field( index{:}, :, 1, 2 ) );
        
        clear dfield d2field
        
        % Get the estimates of the covariances
        VY   = var( Y,  0, D+1 );
        
        CdYY = NaN * ones( [ sY 2] );
        CdYY( index{:}, 1 ) = vectcov( Y, Yx, D+1, 1 );
        CdYY( index{:}, 2 ) = vectcov( Y, Yy, D+1, 1 );
        
        Cd2YY = NaN * ones( [ sY 2 2] );
        Cd2YY( index{:}, 1, 1 ) = vectcov( Y, Yxx, D+1, 1 );
        Cd2YY( index{:}, 1, 2 ) = vectcov( Y, Yxy, D+1, 1 );
        Cd2YY( index{:}, 2, 2 ) = vectcov( Y, Yyy, D+1, 1 );
        Cd2YY( index{:}, 2, 1 ) = Cd2YY( index{:}, 1, 2 );

        CdYdY = NaN * ones( [ sY 3 3] );
        CdYdY( index{:}, 1, 1 ) = var( Yx,  0, D+1 );
        CdYdY( index{:}, 2, 2 ) = var( Yy,  0, D+1 );
        CdYdY( index{:}, 1, 2 ) = vectcov( Yx, Yy, D+1, 1 );
        CdYdY( index{:}, 2, 1 ) = CdYdY( index{:}, 1, 2 );

        Cd2YdY = NaN * ones( [ sY 2 2 2] );
        Cd2YdY( index{:}, 1, 1, 1 ) = vectcov( Yxx, Yx, D+1, 1 );
        Cd2YdY( index{:}, 1, 1, 2 ) = vectcov( Yxx, Yy, D+1, 1 );
        Cd2YdY( index{:}, 2, 2, 1 ) = vectcov( Yyy, Yx, D+1, 1 );
        Cd2YdY( index{:}, 2, 2, 2 ) = vectcov( Yyy, Yy, D+1, 1 );        
        Cd2YdY( index{:}, 1, 2, 1 ) = vectcov( Yxy, Yx, D+1, 1 );
        Cd2YdY( index{:}, 1, 2, 2 ) = vectcov( Yxy, Yy, D+1, 1 );
        Cd2YdY( index{:}, 2, 1, 1 ) = Cd2YdY( index{:}, 1, 2, 1 );
        Cd2YdY( index{:}, 2, 1, 2 ) = Cd2YdY( index{:}, 1, 2, 2 );
        

        % Entries of Christoffel symbol matrix
        for k = 1:D
            for d = 1:D
                for dd = 1:D
                    Gamma.field( index{:}, k, d, dd ) = ...
                          Cd2YdY( index{:}, k, d, dd ) ./ VY ...
                        - CdYY( index{:}, k ) .*...
                                        CdYdY( index{:}, d, dd ) ./ VY.^2 ...
                        - Cd2YY( index{:}, k, d ) .*...
                                        CdYY( index{:}, dd ) ./ VY.^2 ...
                        - CdYY( index{:}, d ) .*...
                                        CdYdY( index{:}, k, dd ) ./ VY.^2 ...
                        + 2 * CdYY( index{:}, k ) .* CdYY( index{:}, d ) .* ...
                                        CdYY( index{:}, dd ) ./  VY.^3;      
                end
            end
        end        
    case 3
        % Rename the partial derivatives of the convolution field
        Yx = squeeze( dfield( index{:}, :, 1 ) );
        Yy = squeeze( dfield( index{:}, :, 2 ) );
        Yz = squeeze( dfield( index{:}, :, 3 ) );
        
        Yxx = squeeze( d2field( index{:}, :, 1, 1 ) );
        Yyy = squeeze( d2field( index{:}, :, 2, 2 ) );
        Yzz = squeeze( d2field( index{:}, :, 3, 3 ) );

        Yxy = squeeze( d2field( index{:}, :, 1, 2 ) );
        Yxz = squeeze( d2field( index{:}, :, 1, 3 ) );
        Yyz = squeeze( d2field( index{:}, :, 2, 3 ) );
        
        clear dfield d2field
        
        % Get the estimates of the covariances
        VY   = var( Y,  0, D+1 );
        
        CdYY = NaN * ones( [ sY 3] );
        CdYY( index{:}, 1 ) = vectcov( Y, Yx, D+1, 1 );
        CdYY( index{:}, 2 ) = vectcov( Y, Yy, D+1, 1 );
        CdYY( index{:}, 3 ) = vectcov( Y, Yz, D+1, 1 );
        
        Cd2YY = NaN * ones( [ sY 3 3] );
        Cd2YY( index{:}, 1, 1 ) = vectcov( Y, Yxx, D+1, 1 );
        Cd2YY( index{:}, 1, 2 ) = vectcov( Y, Yxy, D+1, 1 );
        Cd2YY( index{:}, 1, 3 ) = vectcov( Y, Yxz, D+1, 1 );
        Cd2YY( index{:}, 2, 2 ) = vectcov( Y, Yyy, D+1, 1 );
        Cd2YY( index{:}, 2, 3 ) = vectcov( Y, Yyz, D+1, 1 );
        Cd2YY( index{:}, 3, 3 ) = vectcov( Y, Yzz, D+1, 1 );
        Cd2YY( index{:}, 3, 1 ) = Cd2YY( index{:}, 1, 3 );
        Cd2YY( index{:}, 3, 2 ) = Cd2YY( index{:}, 2, 3 );
        Cd2YY( index{:}, 2, 1 ) = Cd2YY( index{:}, 1, 2 );

        CdYdY = NaN * ones( [ sY 3 3] );
        CdYdY( index{:}, 1, 1 ) = var( Yx,  0, D+1 );
        CdYdY( index{:}, 2, 2 ) = var( Yy,  0, D+1 );
        CdYdY( index{:}, 3, 3 ) = var( Yz,  0, D+1 );
        CdYdY( index{:}, 1, 2 ) = vectcov( Yx, Yy, D+1, 1 );
        CdYdY( index{:}, 1, 3 ) = vectcov( Yx, Yz, D+1, 1 );
        CdYdY( index{:}, 2, 3 ) = vectcov( Yy, Yz, D+1, 1 );
        CdYdY( index{:}, 2, 1 ) = CdYdY( index{:}, 1, 2 );
        CdYdY( index{:}, 3, 1 ) = CdYdY( index{:}, 1, 3 );
        CdYdY( index{:}, 3, 2 ) = CdYdY( index{:}, 2, 3 );

        Cd2YdY = NaN * ones( [ sY 3 3 3] );
        Cd2YdY( index{:}, 1, 1, 1 ) = vectcov( Yxx, Yx, D+1, 1 );
        Cd2YdY( index{:}, 1, 1, 2 ) = vectcov( Yxx, Yy, D+1, 1 );
        Cd2YdY( index{:}, 1, 1, 3 ) = vectcov( Yxx, Yz, D+1, 1 );
        Cd2YdY( index{:}, 2, 2, 1 ) = vectcov( Yyy, Yx, D+1, 1 );
        Cd2YdY( index{:}, 2, 2, 2 ) = vectcov( Yyy, Yy, D+1, 1 );
        Cd2YdY( index{:}, 2, 2, 3 ) = vectcov( Yyy, Yz, D+1, 1 );
        Cd2YdY( index{:}, 3, 3, 1 ) = vectcov( Yzz, Yx, D+1, 1 );
        Cd2YdY( index{:}, 3, 3, 2 ) = vectcov( Yzz, Yy, D+1, 1 );
        Cd2YdY( index{:}, 3, 3, 3 ) = vectcov( Yzz, Yz, D+1, 1 );
        
        Cd2YdY( index{:}, 1, 2, 1 ) = vectcov( Yxy, Yx, D+1, 1 );
        Cd2YdY( index{:}, 1, 2, 2 ) = vectcov( Yxy, Yy, D+1, 1 );
        Cd2YdY( index{:}, 1, 2, 3 ) = vectcov( Yxy, Yz, D+1, 1 );
        Cd2YdY( index{:}, 2, 1, 1 ) = Cd2YdY( index{:}, 1, 2, 1 );
        Cd2YdY( index{:}, 2, 1, 2 ) = Cd2YdY( index{:}, 1, 2, 2 );
        Cd2YdY( index{:}, 2, 1, 3 ) = Cd2YdY( index{:}, 1, 2, 3 );
        
        Cd2YdY( index{:}, 1, 3, 1 ) = vectcov( Yxz, Yx, D+1, 1 );
        Cd2YdY( index{:}, 1, 3, 2 ) = vectcov( Yxz, Yy, D+1, 1 );
        Cd2YdY( index{:}, 1, 3, 3 ) = vectcov( Yxz, Yz, D+1, 1 );
        Cd2YdY( index{:}, 3, 1, 1 ) = Cd2YdY( index{:}, 1, 3, 1 );
        Cd2YdY( index{:}, 3, 1, 2 ) = Cd2YdY( index{:}, 1, 3, 2 );
        Cd2YdY( index{:}, 3, 1, 3 ) = Cd2YdY( index{:}, 1, 3, 3 );
        
        Cd2YdY( index{:}, 2, 3, 1 ) = vectcov( Yyz, Yx, D+1, 1 );
        Cd2YdY( index{:}, 2, 3, 2 ) = vectcov( Yyz, Yy, D+1, 1 );
        Cd2YdY( index{:}, 2, 3, 3 ) = vectcov( Yyz, Yz, D+1, 1 );
        Cd2YdY( index{:}, 3, 2, 1 ) = Cd2YdY( index{:}, 2, 3, 1 );
        Cd2YdY( index{:}, 3, 2, 2 ) = Cd2YdY( index{:}, 2, 3, 2 );
        Cd2YdY( index{:}, 3, 2, 3 ) = Cd2YdY( index{:}, 2, 3, 3 );

        % Entries of Christoffel symbol matrix
        for k = 1:D
            for d = 1:D
                for dd = 1:D
                    Gamma.field( index{:}, k, d, dd ) = ...
                          Cd2YdY( index{:}, k, d, dd ) ./ VY ...
                        - CdYY( index{:}, k ) .*...
                                        CdYdY( index{:}, d, dd ) ./ VY.^2 ...
                        - Cd2YY( index{:}, k, d ) .*...
                                        CdYY( index{:}, dd ) ./ VY.^2 ...
                        - CdYY( index{:}, d ) .*...
                                        CdYdY( index{:}, k, dd ) ./ VY.^2 ...
                        + 2 * CdYY( index{:}, k ) .* CdYY( index{:}, d ) .* ...
                                        CdYY( index{:}, dd ) ./  VY.^3;      
                end
            end
        end        
end

% Remove NaNs due to division by zero
Gamma.field(isnan(Gamma.field(:))) = 0;
%%%%%% END compute Christoffel symbols

return