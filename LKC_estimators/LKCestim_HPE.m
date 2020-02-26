function LKC = LKCestim_HPE( Y, D, mask, Mboot, version )
% LKCestim_HPE( Y, D, mask, Mboot, version )
% This function computes the Lipschitz Killing curvature using the
% Hermite projection estimator proposed in Telschow et al (2020+).
% It uses a fast and exact way to compute the EC curves by looping through
% discrete critical values and computes there contribution to the change in
% Euler characteristic.
%--------------------------------------------------------------------------
% ARGUMENTS
% Y         an array of dimension T_1 x...x T_D x N of N residual
%           fields over an T_1 x...x T_D-square.
%           Note, it is assumed that the fields are standardized residuals
%           or at least mean zero, if Mboot > 1.
% D         dimension of domain of the field
% mask      a boolean array of dimension T_1 x...x T_D. Voxels outside the
%           mask will be set to +oo.
% Mboot     an integer specifying the number of bootstraps used for
%           estimation of LKC. If 1 the estimator is equal to the HPE
%           otherwise the bHPE (recommended Mboot>3e3)
% version   if "C" (default), critical values are computed using C++, if
%           it is any other value a slow matlab only implementation is used.
%--------------------------------------------------------------------------
% OUTPUT
% LKC       the Lipschitz Killing curvatures of the field
%--------------------------------------------------------------------------
% EXAMPLES
% dim   = [50 50 50];
% nsubj = 30;
% D     = length( dim );
% uvals = -6 : 0.1 : 6;
% mask  = true(dim);
% 
% % generate random noise
% Y = noisegen( dim, nsubj, 5*2*sqrt(2*log(2)), 0 );
% % compute HPE and bHPE
% HPE  = LKCestim_HPE( Y, D, mask, 1, "C" );
% bHPE = LKCestim_HPE( Y, D, mask, 5e3, "C" );
% 
% --------------------------------------------------------------------------
% AUTHOR: Fabian Telschow

%%%%%%%%%%%%%%%%%%%%%%% Get constants from the input %%%%%%%%%%%%%%%%%%%%%%
sY     = size( Y );                 % get size of the input array
if length(sY) == D
    sY = [ sY 1 ];
end
N      = sY( D + 1 );               % number of subjects/samples
index  = repmat( {':'}, 1, D );     % get variable domain counter

% check that method is implemented for dimension D
if D > 3
    error( 'Size of Y can only have length<4. Higher dimensional domains not implemented')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% add default values %%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist( 'Mboot', 'var' ) || ~isnumeric(Mboot)
   % default number of bootstrap replicates
   Mboot = 0;
end

if exist( 'mask', 'var' )
    % Check whether mask has the correct size and set values outside the
    % mask to -oo
    if ~all( size( mask ) == sY( 1 : D ) )
        error( 'Incompatible input: The mask needs to have the same size as the first D dimensions of Y.' )
    else
        for i = 1 : N
            tmp = Y( index{:}, i );
            tmp( ~mask ) = -Inf;
            Y( index{:}, i ) = tmp;
        end
    end
else
    % if no mask provided assume the domain is the whole image
    mask = true( sY( 1:D ) );
end
clear i

if nargin < 7
    % default value of "version"
    version = "C";
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the LKC output
if Mboot > 1
    L_hat = zeros( [ D, Mboot ] );
else
    L_hat = zeros( [ D, N ] );    
end

% scaling vector for integral
p = [ sqrt( 2 * pi ); pi; ( 2 * pi )^( 3 / 2 ) / factorial( 3 ); ...
      ( 2 * pi )^( 4 / 2 ) / factorial( 4 ) ];

% Compute LKCs depending on method
if( Mboot > 1 )    
    % Get weights for the multiplier bootstrap
    multiplier = normrnd( 0, 1, [ N, Mboot ] );

    % reshape and and standardize the field, such that it has unit variance
    Y = reshape( Y, prod( sY(1:end-1) ), N );
    % normalize the residuals
    Y = Y ./ sqrt(sum( Y.^2, 2 ));

    for i = 1:Mboot
        % get the bootstrapped process
        mY = reshape( Y * multiplier( :, i ), sY(1:end-1) );

        % Get the EC stepfunctions
        EC = EulerCharCrit( mY, D, mask, version );
        EC = EC{ 1 };

        % Get LKC by integrating the EC curves against the Hermite
        % polynomials
        v =  EC( 2:end-1, 1 )';
        % differences of EC between consecutive critical values
        b = -diff( EC( 2:end, 2 ) );
        % Hermite polynomials evaluated at critical values
        H =  [ v; ( v.^2 - 1 ); ( v.^3 - 3*v ); ( v.^4 - 6*v.^2 + 3 ) ];

        L_hat( :, i ) = p( 1:D ) .* ( H( 1:D, : ) * b );
    end
    L0 = EC(1,2);
else
    % Get the EC stepfunctions
    ECall = EulerCharCrit( Y, D, mask, version );

    for i = 1:N
        % Get LKC by integrating the Euler Char curves against the Hermite
        % polynomials
        v =  ECall{ i }( 2:end-1, 1 )';
        % differences of EC between consecutive critical values
        b = -diff( ECall{ i }( 2:end, 2 ) );
        % Hermite polynomials evaluated at critical values
        H =  [ v; ( v.^2 - 1 ); ( v.^3 - 3 * v ); ( v.^4 - 6 * v.^2 + 3 ) ];

        L_hat( :, i ) = p( 1:D ) .* ( H( 1:D, : ) * b );
    end
    L0 = ECall{1}(1,2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% stat summary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Summary estimators: LKCs
if N > 1
    L_hatn     = mean( L_hat, 2 );
    Sigma_hat  = cov( L_hat' );
    L_se_hat   = sqrt( diag( Sigma_hat ) / size( L_hat, 2 ) );
    L_conf_hat = cat( 2, L_hatn - 1.96 * L_se_hat, ...
                         L_hatn + 1.96 * L_se_hat );
else
    L_hatn     = L_hat;
    Sigma_hat  = NaN;
    L_se_hat   = NaN;
    L_conf_hat = NaN;
end
% Summarize output
LKC  = struct( 'hat1', L_hat, 'hatn', L_hatn, 'Sigma_hat', Sigma_hat, ...
               'se_hat', L_se_hat, 'conf_hat', L_conf_hat, 'L0', L0 );
return