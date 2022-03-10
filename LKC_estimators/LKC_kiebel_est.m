function [ L, Lambda, FWHM, L_forman, FWHM_forman ] = LKC_kiebel_est( field, k )
% LKC_kiebel_est( field, k ) estimates the LKCs assuming isotropy according
% to the article [1].
%
% [1] Kiebel, Stefan J., et al. "Robust smoothness estimation in
% statistical parametric maps using standardized residuals from the
% general linear model." Neuroimage 10.6 (1999): 756-766.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  field   an object of class Field representing observations of a field,
%          fiberD = 1 and fibersize > 1.
% Optional
%  k      integer such that df = field.fieldsize - k. It is the total
%         degree of freedoms of the input Field object, which in general
%         should be residuals, e.g., from a linear model, compare Kiebel
%         1999 eq. (10).
%--------------------------------------------------------------------------
% OUTPUT
%  L   an 1 x field.D vector containing the LKCs L_1,...,L_field.D
%  Lambda  a D x D matrix of the estimated Riemannian metric. Note this is
%          assumed constant across the domain, since Kiebel assumes
%          stationarity.
%  FWHM    a 1 x D vector of estimates of FWHM based on the Lambda matrix 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow, Samuel Davenport
%--------------------------------------------------------------------------

%% Check mandatory input
%--------------------------------------------------------------------------
% Get the important variables from the Field object
Y      = field.field;
mask   = zero2nan( field.mask );
D      = field.D;
voxdim = get_dx( field );
nsubj  = field.fibersize;

index = repmat( {':'}, [ 1 D ] );

% check that method is implemented for dimension D
if( D > 4 )
    error('The method is currently only implemented for field domains of dimension D<5.')
end

%% Main function
%--------------------------------------------------------------------------
% Get the Standardized residuals eq (13)
normY = sqrt( sum( Y.^2, D + 1) ); 

% Mask the residuals
Y = ( Y ./ normY ) .* mask;

% Factor from eq. (14) Kienel 1999
dffac = ( nsubj - k - 2 ) / ( nsubj - k - 1 );

% Preallocate the Lambda matrix
Lambda = zeros( D );

switch D
    case 1
        % Symmetric differential quotient eq. (15) Kiebel
        Xderiv = ( Y(3:end,:) - Y(1:end-2,:) ) / voxdim(1) / 2;

        % Get the variance of this derivative
        tmp    = ~isnan( Xderiv( index{:}, 1 ) );
        Lambda(1,1) = sum( Xderiv( ~isnan( Xderiv ) ).^2 )...
                                                / sum( tmp(:) ) * dffac;
    case 2
        % Symmetric differential quotient eq. (15) Kiebel first component
        Xderiv = ( Y(3:end,2:end-1,:) - Y(1:end-2,2:end-1,:) )...
                                                    / voxdim(1) / 2;

        % Get the variance of this derivative
        tmp    = ~isnan( Xderiv( index{:}, 1 ) );
        Lambda(1,1) = sum( Xderiv( ~isnan( Xderiv ) ).^2 )...
                                                / sum( tmp(:) ) * dffac;
                                            
        % Symmetric differential quotient yeq. (15) Kiebel second component
        Yderiv = ( Y(2:end-1,3:end,:) - Y(2:end-1,1:end-2,:) )...
                                                    / voxdim(2) / 2;

        % Get the variance of this derivative
        tmp    = ~isnan( Yderiv( index{:}, 1 ) );
        Lambda(2,2) = sum( Yderiv( ~isnan( Yderiv ) ).^2 )...
                                                / sum( tmp(:) ) * dffac;
                                            
        % Cross term
        Yderiv = Xderiv.*Yderiv;
        tmp     = ~isnan( Yderiv( index{:}, 1 ) );
        Lambda(1,2) = sum( Yderiv( ~isnan( Yderiv ) ) )...
                                                / sum( tmp(:) ) * dffac;
        Lambda(1,2) = Lambda(2,1);

    case 3
        % Symmetric differential quotient eq. (15) Kiebel first component
        Xderiv = ( Y(3:end,2:end-1,2:end-1,:) - Y(1:end-2,2:end-1,2:end-1,:) )...
                                                    / voxdim(1) / 2;

        % Get the variance of this derivative
        tmp    = ~isnan( Xderiv( index{:}, 1 ) );
        Lambda(1,1) = sum( Xderiv( ~isnan( Xderiv ) ).^2 )...
                                                / sum( tmp(:) ) * dffac;
                                            
        % Symmetric differential quotient yeq. (15) Kiebel second component
        Yderiv = ( Y(2:end-1,3:end,2:end-1,:) - Y(2:end-1,1:end-2,2:end-1,:) )...
                                                    / voxdim(2) / 2;

        % Get the variance of this derivative
        tmp    = ~isnan( Yderiv( index{:}, 1 ) );
        Lambda(2,2) = sum( Yderiv( ~isnan( Yderiv ) ).^2 )...
                                                / sum( tmp(:) ) * dffac;
                                            
        % Cross term
        XYderiv = Xderiv.*Yderiv;
        tmp     = ~isnan( XYderiv( index{:}, 1 ) );
        Lambda(1,2) = sum( XYderiv( ~isnan( XYderiv ) ) )...
                                                / sum( tmp(:) ) * dffac;
        Lambda(2,1) = Lambda(1,2);
        
        clear XYderiv
        
        % Symmetric differential quotient yeq. (15) Kiebel third component
        Zderiv = ( Y(2:end-1,2:end-1,3:end,:) - Y(2:end-1,2:end-1,1:end-2,:) )...
                                                    / voxdim(3) / 2;

        % Get the variance of this derivative
        tmp    = ~isnan( Zderiv( index{:}, 1 ) );
        Lambda(3,3) = sum( Zderiv( ~isnan( Zderiv ) ).^2 )...
                                                / sum( tmp(:) ) * dffac;
                                            
        % Cross termS
        XZderiv = Xderiv.*Zderiv;
        tmp     = ~isnan( XZderiv( index{:}, 1 ) );
        Lambda(1,3) = sum( XZderiv( ~isnan( XZderiv ) ) )...
                                                / sum( tmp(:) ) * dffac;
        Lambda(3,1) = Lambda(1,3);
        
        % Cross termS
        YZderiv = Yderiv.*Zderiv;
        tmp     = ~isnan( YZderiv( index{:}, 1 ) );
        Lambda(2,3) = sum( YZderiv( ~isnan( YZderiv ) ) )...
                                                / sum( tmp(:) ) * dffac;
        Lambda(3,2) = Lambda(2,3);
        
end

FWHM =  sqrt( 4*log(2) ./ diag(Lambda) );
voxmfd  = VoxManifold( constfield( Lambda, field.mask, field.xvals ) );
% Obtain the LKCs
L = LKC_est( voxmfd );

%--------------------------------------------------------------------------
% Forman estimators
sigma_est_forman = sqrt( -1 ./ 4 ./ log( 1 - diag(Lambda) / 2 ) );
FWHM_forman  = sigma2FWHM( sigma_est_forman );

Lambda_f = 4*log(2) * diag( 1 ./ FWHM_forman.^2 );

% Obtain the LKCs from Forman FWHM
voxmfd    = VoxManifold( constfield( Lambda_f, field.mask, field.xvals ) );
L_forman  = LKC_est( voxmfd );
%--------------------------------------------------------------------------



