function Lambda = Lambda_lat_stat_est( field, num_deriv, nu, pool )
% Lambda_lat_stat_est( field, num_deriv, scale ) estimates the Riemmanian
% metric under the  assumption of stationarity from the lattice values.
% This is essentially the approach taken in [1] and [2]. 
%
% [1] Kiebel, Stefan J., et al. "Robust smoothness estimation in
% statistical parametric maps using standardized residuals from the
% general linear model." Neuroimage 10.6 (1999): 756-766.
%
% [2] Forman, Steven D., et al. "Improved assessment of significant
% activation in functional magnetic resonance imaging (fMRI): Use of a
% clusterâ€?size threshold." Magnetic Resonance in medicine 33.5 (1995): 636-647.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  field      an object of class Field representing observations of a field,
%             fiberD = 1 and fibersize > 1.
% Optional
%  num_deriv  a string. Choices are "symmetric" or "onesided" for symmetric
%             or onesided numeric derivatives. Default is "symmetric".
%  scale      a numeric. If scale ~= 0 the estimate is multiplied by
%             ( nsubj - scale - 2 ) / ( nsubj - scale - 1 ).
%             If scale = 0, then no constant is multiplied to the estimate.
%             Default is 1.
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Fill the missing variables
if ~exist('num_deriv', 'var')
   num_deriv = "symmetric"; 
end

if ~exist('nu', 'var')
    nsubj  = field.fibersize;
    % Kiebel 99: eq. 10
    nu = nsubj;
end

if ~exist('pool', 'var')
   pool = false; 
end

% Check the input of num_deriv
if ~(strcmp(num_deriv, "symmetric") || strcmp(num_deriv, "onesided"))
   error("num_deriv must be either 'onesided' or 'symmetric'");
end


% Get the important variables from the Field object
Y      = field.field;
mask   = zero2nan( field.mask );
D      = field.D;
voxdim = get_dx( field );

index = repmat( {':'}, [ 1 D ] );

% check that method is implemented for dimension D
if( D > 4 )
    error('The method is currently only implemented for field domains of dimension D<5.')
end


%% Main function
%--------------------------------------------------------------------------
if ~pool
    % Get residuals eq. (6)
    Y = Y - mean( Y, D+1 );

    % Get the norm of the Standardized residuals: denominator eq (13)
    normY = sqrt( sum( Y.^2, D + 1) ); 

    % Compute and Mask the standardized residuals
    Y = ( Y ./ normY ) .* mask;

    % Factor from eq. (14) Kiebel 1999, N computed later as sum(tmp(:))!
    dffac = (nu - 2) / (nu - 1);
    dffac = 1;
else
    Y = Y - mean( Y, D+1 );
    % Estimate the variance across the image
    varY = var(field);
    varY = mean(varY.field(varY.field(:)~=0));
    
    % Compute and Mask the standardized residuals
    Y = Y / sqrt(varY) .* mask;
    
    % Factor for unbiased estimation of the variance of the derivatives
    % we assume that varY is estimated precisely due to averaging over the
    % whole domain
    dffac = 1 / (field.fibersize - 1);
end

% Preallocate the Lambda matrix
Lambda = zeros( D );

switch num_deriv
    case "symmetric"
        switch D
            case 1
                % Symmetric differential quotient eq. (15) Kiebel
                Xderiv = ( Y(3:end,:) - Y(1:end-2,:) ) / voxdim(1) / 2;

                % Get the variance of this derivative
                tmp    = ~isnan( Xderiv( index{:}, 1 ) );
                Lambda(1,1) = sum( Xderiv( ~isnan( Xderiv ) ).^2 )...
                                                / sum( tmp(:) );
            case 2
                % Symmetric differential quotient eq. (15) Kiebel first component
                Xderiv = ( Y(3:end,:,:) - Y(1:end-2,:,:) )...
                                                    / voxdim(1) / 2;

                % Get the variance of this derivative
                tmp    = ~isnan( Xderiv( index{:}, 1 ) );
                Lambda(1,1) = sum( Xderiv( ~isnan( Xderiv ) ).^2 )...
                                                / sum( tmp(:) );
                                            
                % Symmetric differential quotient eq. (15) Kiebel second component
                Yderiv = ( Y(:,3:end,:) - Y(:,1:end-2,:) )...
                                                    / voxdim(2) / 2;

                % Get the variance of this derivative
                tmp    = ~isnan( Yderiv( index{:}, 1 ) );
                Lambda(2,2) = sum( Yderiv( ~isnan( Yderiv ) ).^2 )...
                                                / sum( tmp(:) );
                                            
                % Cross term
                Yderiv = Xderiv(:, 2:end-1, :) .* Yderiv(2:end-1, :, :);
                tmp     = ~isnan( Yderiv( index{:}, 1 ) );
                Lambda(1,2) = sum( Yderiv( ~isnan( Yderiv ) ) )...
                                                / sum( tmp(:) );
                Lambda(1,2) = Lambda(2,1);

            case 3
                % Symmetric differential quotient eq. (15) Kiebel first component
                Xderiv = ( Y(3:end,:,:,:) - Y(1:end-2,:,:,:) )...
                                                    / voxdim(1) / 2;

                % Get the variance of this derivative
                tmp    = ~isnan( Xderiv( index{:}, 1 ) );
                Lambda(1,1) = sum( Xderiv( ~isnan( Xderiv ) ).^2 )...
                                                / sum( tmp(:) );
                                            
                % Symmetric differential quotient yeq. (15) Kiebel second component
                Yderiv = ( Y(:, 3:end, :, :) - Y(:, 1:end-2, :, :) )...
                                                    / voxdim(2) / 2;

                % Get the variance of this derivative
                tmp    = ~isnan( Yderiv( index{:}, 1 ) );
                Lambda(2,2) = sum( Yderiv( ~isnan( Yderiv ) ).^2 )...
                                                / sum( tmp(:) );
                                            
                % Cross term
                XYderiv = Xderiv(:, 2:end-1, 2:end-1, :) .* Yderiv(2:end-1, :, 2:end-1,:);
                tmp     = ~isnan( XYderiv( index{:}, 1 ) );
                Lambda(1,2) = sum( XYderiv( ~isnan( XYderiv ) ) )...
                                                / sum( tmp(:) );
                Lambda(2,1) = Lambda(1,2);
        
                clear XYderiv
        
                % Symmetric differential quotient yeq. (15) Kiebel third component
                Zderiv = ( Y(:, :, 3:end, :) - Y(:, :, 1:end-2, :) )...
                                                    / voxdim(3) / 2;

                % Get the variance of this derivative
                tmp    = ~isnan( Zderiv( index{:}, 1 ) );
                Lambda(3,3) = sum( Zderiv( ~isnan( Zderiv ) ).^2 )...
                                                / sum( tmp(:) );
                                            
                % Cross termS
                XZderiv = Xderiv(:, 2:end-1, 2:end-1, :) .* Zderiv(2:end-1, 2:end-1, :, :);
                tmp     = ~isnan( XZderiv( index{:}, 1 ) );
                Lambda(1,3) = sum( XZderiv( ~isnan( XZderiv ) ) )...
                                                / sum( tmp(:) );
                Lambda(3,1) = Lambda(1,3);
        
                % Cross termS
                YZderiv = Yderiv(2:end-1, :, 2:end-1,:) .* Zderiv(2:end-1, 2:end-1, :, :);
                tmp     = ~isnan( YZderiv( index{:}, 1 ) );
                Lambda(2,3) = sum( YZderiv( ~isnan( YZderiv ) ) )...
                                                / sum( tmp(:) );
                Lambda(3,2) = Lambda(2,3);
        
        end
        
    case "onesided"
        switch D
            case 1
                % Symmetric differential quotient eq. (15) Kiebel
                Xderiv = ( Y(2:end,:) - Y(1:end-1,:) ) / voxdim(1);

                % Get the variance of this derivative
                tmp    = ~isnan( Xderiv( index{:}, 1 ) );
                Lambda(1,1) = sum( Xderiv( ~isnan( Xderiv ) ).^2 )...
                                                / sum( tmp(:) );
            case 2
                % Symmetric differential quotient eq. (15) Kiebel first component
                Xderiv = ( Y(2:end,:,:) - Y(1:end-1,:,:) )...
                                                    / voxdim(1);

                % Get the variance of this derivative
                tmp    = ~isnan( Xderiv( index{:}, 1 ) );
                Lambda(1,1) = sum( Xderiv( ~isnan( Xderiv ) ).^2 )...
                                                / sum( tmp(:) );
                                            
                % Symmetric differential quotient eq. (15) Kiebel second component
                Yderiv = ( Y(:,2:end,:) - Y(:,1:end-1,:) )...
                                                    / voxdim(2);

                % Get the variance of this derivative
                tmp    = ~isnan( Yderiv( index{:}, 1 ) );
                Lambda(2,2) = sum( Yderiv( ~isnan( Yderiv ) ).^2 )...
                                                / sum( tmp(:) );
                                            
                % Cross term
                Yderiv = Xderiv(:,1:end-1,:) .* Yderiv(1:end-1, :, :);
                tmp     = ~isnan( Yderiv( index{:}, 1 ) );
                Lambda(1,2) = sum( Yderiv( ~isnan( Yderiv ) ) )...
                                                / sum( tmp(:) );
                Lambda(1,2) = Lambda(2,1);

            case 3
                % Symmetric differential quotient eq. (15) Kiebel first component
                Xderiv = ( Y(2:end,:,:,:) - Y(1:end-1,:,:,:) )...
                                                    / voxdim(1);

                % Get the variance of this derivative
                tmp    = ~isnan( Xderiv( index{:}, 1 ) );
                Lambda(1,1) = sum( Xderiv( ~isnan( Xderiv ) ).^2 )...
                                                / sum( tmp(:) );
                                            
                % Symmetric differential quotient yeq. (15) Kiebel second component
                Yderiv = ( Y(:,2:end,:,:) - Y(:,1:end-1,:,:) )...
                                                    / voxdim(2);

                % Get the variance of this derivative
                tmp    = ~isnan( Yderiv( index{:}, 1 ) );
                Lambda(2,2) = sum( Yderiv( ~isnan( Yderiv ) ).^2 )...
                                                / sum( tmp(:) );
                                            
                % Cross term
                XYderiv = Xderiv(:,1:end-1,1:end-1,:).*Yderiv(1:end-1,:,1:end-1,:);
                tmp     = ~isnan( XYderiv( index{:}, 1 ) );
                Lambda(1,2) = sum( XYderiv( ~isnan( XYderiv ) ) )...
                                                / sum( tmp(:) );
                Lambda(2,1) = Lambda(1,2);
        
                clear XYderiv
        
                % Symmetric differential quotient yeq. (15) Kiebel third component
                Zderiv = ( Y(:,:,2:end,:) - Y(:,:,1:end-1,:) )...
                                                    / voxdim(3);

                % Get the variance of this derivative
                tmp    = ~isnan( Zderiv( index{:}, 1 ) );
                Lambda(3,3) = sum( Zderiv( ~isnan( Zderiv ) ).^2 )...
                                                / sum( tmp(:) );
                                            
                % Cross termS
                XZderiv = Xderiv(:, 1:end-1, 1:end-1, :) .* Zderiv(1:end-1, 1:end-1, :, :);
                tmp     = ~isnan( XZderiv( index{:}, 1 ) );
                Lambda(1,3) = sum( XZderiv( ~isnan( XZderiv ) ) )...
                                                / sum( tmp(:) );
                Lambda(3,1) = Lambda(1,3);
        
                % Cross termS
                YZderiv = Yderiv(1:end-1, :, 1:end-1, :) .* Zderiv(1:end-1, 1:end-1, :, :);
                tmp     = ~isnan( YZderiv( index{:}, 1 ) );
                Lambda(2,3) = sum( YZderiv( ~isnan( YZderiv ) ) )...
                                                / sum( tmp(:) );
                Lambda(3,2) = Lambda(2,3);
        
        end
end

Lambda = Lambda * dffac;

return
