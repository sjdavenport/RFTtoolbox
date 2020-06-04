function LKC = LKC_HP_est( Y, mask, Mboot, normalize, version )
% LKC_HPE_est( Y, mask, Mboot, normalize, version ) computes the Lipschitz
% Killing curvature using the Hermite projection estimator proposed in
% Telschow et al (2020+).
% It uses a fast and exact way to compute the EC curves by looping through
% discrete critical values and computes their contribution to the change in
% the Euler characteristic.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory:
%  Y      an array of dimension T_1 x...x T_D x N of N random
%         fields over an T_1 x...x T_D rectangle.
%  mask   a logical array of dimension T_1 x...x T_D.
%
% Optional:
%  Mboot     an integer specifying the number of bootstraps used for
%            estimation of LKC. If "1" the HPE is used otherwise the
%            bHPE. Default 1.
%  normalize logical indicating whether Y needs to be standardized. 
%            Default 1, i.e., mean will be subtracted and data will be 
%            standardized to have empirical variance 1, if N>1 else 0.
%  version   if "C" (default), Euler characteristic is computed using a
%            routine to compute critical values in C++, else a slow matlab
%            only implementation is used, which is not recommended for D>2.
%--------------------------------------------------------------------------
% OUTPUT
%  LKC     structure containing fields:
%          - hatL1: D x N or D x Mboot array containing the LKC estimates
%                   for each random field
%          - hatL: D x 1 vector of estimates of LKC for the sample Y. It
%                  is the average of hatL1.
%          - L0: integer containing the Euler characteristic of the mask
%                equivalently the zeroth LKC of the random fields.
%          - hatSIGMA: D x D estimate of the covariance matrix of the
%                      estimate hatL. It is computed using the empirical
%                      covariance of hatL1.
%          - hatSE:    D x 1 vector of standard errors for components of
%                      hatL based on the CLT. 
%          - confInt95: approximate 95% confidence intervals for hatL
%                       based on the standard CLT.
%--------------------------------------------------------------------------
% EXAMPLES
% %%% D = 1
% dim   = [ 100 ];
% mask  = ones( [dim 1] );
% FWHM  = 12;
% 
% % generate random noise
% Y = noisegen( dim, 10, FWHM, 0 );
% 
% % compute HPE for one subject. Note that normalize needs to be 0! It is
% % automatically correctly chosen if not specified. All values give the same
% % result!
% HPE_C = LKCestim_HPE( Y(:,1), mask, 1, 0, "C" );
% HPE_m = LKCestim_HPE( Y(:,1), mask, 1, 0, "matlab" );
% HPE_auto = LKCestim_HPE( Y(:,1), mask );
% 
% % compute HPE for several subject. Here we can choose to normalize, which
% % is recommended.
% HPE_C = LKCestim_HPE( Y, mask, 1, 0, "C" );
% HPE_m = LKCestim_HPE( Y, mask, 1, 1, "matlab" );
% HPE_auto = LKCestim_HPE( Y, mask ); % same as HPE_m using "C"
% 
% % compute HPE for several subject. Here normalize does not have any
% % since it is required to normalize and done automatically. Note that bHPE
% % requires several subjects
% bHPE_C = LKCestim_HPE( Y, mask, 1e3, 0, "C" );
% bHPE_m = LKCestim_HPE( Y, mask, 1e3, 1, "matlab" );
% 
% %%%% D = 2
% dim   = [ 50 50 ];
% mask  = ones( dim );
% FWHM  = sigma2FWHM(5);
% mask(25,25) = 1;
% 
% % generate random noise
% Y = noisegen( dim, 10, FWHM );
% 
% % compute HPE for one subject. Note that normalize needs to be 0! It is
% % automatically correctly chosen if not specified. All values give the same
% % result!
% HPE_C = LKCestim_HPE( Y(:,:,1), mask, 1, 0, "C" );
% HPE_m = LKCestim_HPE( Y(:,:,1), mask, 1, 0, "matlab" );
% HPE_auto = LKCestim_HPE( Y(:,:,1), mask );
% 
% % compute HPE for several subject. Here we can choose to normalize, which
% % is recommended, since usually it is unknown whether the observations
% % have mean zero and unit variance.
% HPE_C = LKCestim_HPE( Y, mask, 1, 0, "C" );
% HPE_m = LKCestim_HPE( Y, mask, 1, 1, "matlab" );
% HPE_auto = LKCestim_HPE( Y, mask ); % same as HPE_m using "C"
% 
% % compute HPE for several subject. Here normalize does not have any
% % since it is required to normalize and done automatically. Note that bHPE
% % requires several subjects. Note value is unexpectedly low. Need to figure
% % out why. Using "matlab" not recommended since it is extremely slow.
% bHPE_C = LKCestim_HPE( Y, mask, 3e3, 0, "C" );
% bHPE_m = LKCestim_HPE( Y, mask, 1e3, 0, "matlab" );
% 
% %%% D = 3
% dim   = [ 35 35 35 ];
% mask  = ones( dim );
% FWHM  = 5;
% mask(25,25,25) = 0;
% 
% % generate random noise
% Y = noisegen( dim, 10, FWHM );
% 
% % compute HPE for one subject. Note that normalize needs to be 0! It is
% % automatically correctly chosen if not specified. All values give the same
% % result!
% HPE_C = LKCestim_HPE( Y(:,:,:,1), mask, 1, 0, "C" );
% HPE_m = LKCestim_HPE( Y(:,:,:,1), mask, 1, 0, "matlab" );
% HPE_auto = LKCestim_HPE( Y(:,:,:,1), mask );
% 
% % compute HPE for several subject. Here we can choose to normalize, which
% % is recommended.
% HPE_C = LKCestim_HPE( Y, mask, 1, 0, "C" );
% HPE_m = LKCestim_HPE( Y, mask, 1, 1, "matlab" );
% HPE_auto = LKCestim_HPE( Y, mask ); % same as HPE_m using "C"
% 
% % compute HPE for several subject. Here normalize does not have any
% % since it is required to normalize and done automatically. Note that bHPE
% % requires several subjects.
% % Never use "matlab" for bHPE since it is extremely slow.
% bHPE_C = LKCestim_HPE( Y, mask, 3e3, 0, "C" );
% -------------------------------------------------------------------------
% AUTHOR: Fabian Telschow
%--------------------------------------------------------------------------


%% %-----------------------------------------------------------------------
%  check mandatory input and get important constants
%--------------------------------------------------------------------------
%%% get constants from the input
% size of the domain
sM = size( mask );
sY = size( Y );

% dimension of the domain, since matlab is not consistent for D<1, we need
% to catch this case 
if sM( 2 ) == 1 && length( sM ) == 2
    D = 1;
else
    D = length(sM);
end

% modify the size array of Y to keep track of N.
if length( sY ) == D
    sY = [ sY 1 ];
end
% get number of subjects/samples
N      = sY( D + 1 );

% check validity of mask input
if ~all( sM == sY( 1:end-1 ) ) && ~all( sM == sY ) && ...
   ~( D == 1 && sM( 1 ) == sY( 1 ) )
   error( 'The mask needs to have the same size as Y. Note that realisations are stored as columns.' )
end

% check that method is implemented for dimension D
if D > 3
    error( 'D must be < 4. Higher dimensional domains have not been implemented')
end


%% %-----------------------------------------------------------------------
%  add and check optional values
%--------------------------------------------------------------------------
if ~exist( 'Mboot', 'var' )
   % default number of bootstrap replicates
   Mboot = 1;
end

if N == 1 && Mboot > 1
    error('The bHPE requires that N > 2!')
end

if ~exist( 'normalize', 'var' )
    % default value of "normalize"
    if N == 1
       normalize = 0;
    else
       normalize = 1;
    end
end

if ~exist( 'version', 'var' )
    % default value of "version"
    version = "C";
end


%% %-----------------------------------------------------------------------
%  main function
%--------------------------------------------------------------------------
% initialize the LKC output
if Mboot > 1
    L_hat = NaN*zeros( [ D, Mboot ] );
else
    L_hat = NaN*zeros( [ D, N ] );    
end

% apply the mask to the data, i.e., set values outside the mask to -oo
if ~all( mask(:) )
    mask( mask == 0 ) = -Inf;
    Y = repmat( mask, [ ones( [ 1 D ] ), N ] ).* Y;
end

% scaling vector for integral
p = [ sqrt( 2 * pi ); pi; ( 2 * pi )^( 3 / 2 ) / factorial( 3 ); ...
      ( 2 * pi )^( 4 / 2 ) / factorial( 4 ) ];

% compute LKCs depending on method
if( Mboot > 1 )    
    % get weights for the multiplier bootstrap
    multiplier = normrnd( 0, 1, [ N, Mboot ] );

    % reshape and and standardize the field, such that it has unit variance
    Y = reshape( Y, prod( sY(1:end-1) ), N );
    % normalize the residuals
    Y = Y - mean( Y, 2 );
    Y = Y ./ sqrt( sum( Y.^2, 2 ) );

    for i = 1:Mboot
        % get the bootstrapped process
        if D~=1
            mY = reshape( Y * multiplier( :, i ), sY(1:end-1) );
        else
            mY = Y * multiplier( :, i );
        end

        % get the EC stepfunctions
        EC = EulerCharCrit( mY, D, mask, version );
        EC = EC{ 1 };

        % get LKC by integrating the EC curves against the Hermite
        % polynomials
        v =  EC( 2:end-1, 1 )';
        % differences of EC between consecutive critical values
        b = -diff( EC( 2:end, 2 ) );
        % Hermite polynomials evaluated at critical values
        H =  [ v; ( v.^2 - 1 ); ( v.^3 - 3*v ); ( v.^4 - 6*v.^2 + 3 ) ];

        L_hat( :, i ) = p( 1:D ) .* ( H( 1:D, : ) * b );
    end
    L0 = EC( 1, 2 );
else
    % normalize the field to have mean zero and unit variance
    if normalize
        Y = ( Y - mean( Y, D+1 ) ) ./ std( Y, 0, D+1 );
    end
    
    % get the EC stepfunctions
    ECall = EulerCharCrit( Y, D, mask, version );

    for i = 1:N
        % get LKC by integrating the Euler Char curves against the Hermite
        % polynomials
        v =  ECall{ i }( 2:end-1, 1 )';
        % differences of EC between consecutive critical values
        b = -diff( ECall{ i }( 2:end, 2 ) );
        % Hermite polynomials evaluated at critical values
        H =  [ v; ( v.^2 - 1 ); ( v.^3 - 3 * v );...
                    ( v.^4 - 6 * v.^2 + 3 ) ];

        L_hat( :, i ) = p( 1:D ) .* ( H( 1:D, : ) * b );
    end
    L0 = ECall{1}(1,2);
end


%% %-----------------------------------------------------------------------
%  prepare output structure
%--------------------------------------------------------------------------
% summary estimators: LKCs as mean of the individual estimators and the
% corresponding covariance matrix and confidence intervals
if N > 1
    L_hatn     = mean( L_hat, 2 );
    Sigma_hat  = cov( L_hat' );
    L_se_hat   = sqrt( diag( Sigma_hat ) / size( L_hat, 2 ) );
    L_conf_hat = cat( 2, L_hatn - 1.96 * L_se_hat, ...
                         L_hatn + 1.96 * L_se_hat );
else
    L_hatn     = L_hat;
    Sigma_hat  = "Covariance matrix cannot be computed for a single field";
    L_se_hat   = "standard error cannot be computed for a single field";
    L_conf_hat = "Confidence intervalls can not be computed for a single field";
end
% summarize output
LKC  = struct(  'hatL1', L_hat, 'hatL', L_hatn, 'L0', L0, 'hatSIGMA',...
                Sigma_hat, 'hatSE', L_se_hat, 'confInt95', L_conf_hat );
            
return