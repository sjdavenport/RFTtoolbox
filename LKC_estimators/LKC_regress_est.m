function [LKC_est, LKC_est_var] = LKC_regress_est( field, uvals, CovEst,...
                                                   L, standardize, Mboot, version )
% LKC_REGRESS_EST( Y, mask, Mboot, normalize, version ) computes the Lipschitz
% Killing curvature using the Hermite projection estimator proposed in
% Telschow et al (2020+).
% It uses a fast and exact way to compute the EC curves by looping through
% discrete critical values and computes their contribution to the change in
% the Euler characteristic.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory:
%  field  an object of class field containing observation of a random field.
%
% Optional:
%  uvals     an vector indicating at which values the EC is evaluated.
%            Default value -3:0.2:3.
%  CovEst    a string giving the method to estimate the covariance.
%            Options are:
%               - OLS: covariance matrix is identity
%  L         a numeric vector of length field. D + 1 indicating the values
%            of known LKCs. NaN indicate that the LKC needs to be
%            estimated. Default [1, NaN(1, D)]
%  standardize logical indicating whether Y needs to be standardized. 
%              Default 1, i.e., mean will be subtracted and data will be 
%              standardized to have empirical variance 1, if N>1 else 0.
%  Mboot     an integer specifying the number of bootstraps used for
%            estimation of LKC. If "1" the HPE is used otherwise the
%            bHPE. Default 1.
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
% %% 1D Isotropic field
% %--------------------------------------------------------------------------
% % Parameters
% D     = 1;
% FWHM  = 6;
% nsubj = 50;
% T     = 100;
% Msim  = 100;
% 
% % Generate mask
% pad    = ceil( 4*FWHM2sigma( FWHM ) );
% mask = pad_vals( true( [ T, 1 ] ), pad, false );
% lat_masked = false;
% 
% % Get high resolution theory
% resThy = 301;
% paramsThy = ConvFieldParams( FWHM, resThy, ceil(resThy/2), lat_masked );
% mask = pad_vals( true( [ T, 1 ] ), pad, false );
% theoryLKC = LKC_wncfield_theory( mask, paramsThy );
% 
% % Generate data
% lat_data =  wfield( mask, nsubj );
% 
% % Convolution fields
% cfield  = convfield( lat_data, params );
% dcfield = convfield( lat_data, params, 1 );
% 
% 
% % Regression estimators
% % values for regression
% uvals = -4:0.5:4;
% 
% % Which LKC needs to be estimated?
% L = [1, NaN];
% 
% L_reg_OLS  = LKC_regress_est( cfield, uvals, "OLS", L, 1, 1 );
% L_breg_OLS = LKC_regress_est( cfield, uvals, "OLS", L, 1, 200 );
% L_reg_SD   = LKC_regress_est( cfield, uvals, "SD", L, 1, 1 );
% L_breg_SD  = LKC_regress_est( cfield, uvals, "SD", L, 1, 200 );
% L_reg_SC   = LKC_regress_est( cfield, uvals, "SC", L, 1, 1 );
% L_breg_SC  = LKC_regress_est( cfield, uvals, "SC", L, 1, 200 );
% L_reg_PI   = LKC_regress_est( cfield, uvals, "PI", L, 1, 1 );
% L_breg_PI  = LKC_regress_est( cfield, uvals, "PI", L, 1, 200 );
% 
% struct(...
%     'theoryL', theoryLKC,...
%     'L_reg_OLS', L_reg_OLS,...
%     'L_breg_OLS', L_reg_OLS,...
%     'L_reg_SD', L_reg_SD,...
%     'L_breg_SD', L_breg_SD,...
%     'L_reg_SC', L_reg_SC,...
%     'L_breg_SC', L_breg_SC,...
%     'L_reg_PI', L_reg_PI,...
%     'L_breg_PI', L_breg_PI...
% )
% 
% %% 2D Isotropic field
% %--------------------------------------------------------------------------
% % Parameters
% D     = 2;
% FWHM  = 6;
% nsubj = 50;
% T     = 100;
% Msim  = 100;
% 
% % Generate mask
% pad    = ceil( 4*FWHM2sigma( FWHM ) );
% mask = pad_vals( true( [ T, T ] ), pad, false );
% lat_masked = false;
% 
% % Get high resolution theory
% resThy = 7;
% paramsThy = ConvFieldParams( [FWHM FWHM], resThy, ceil(resThy/2), lat_masked );
% theoryLKC = LKC_wncfield_theory( mask, paramsThy );
%  
% % Generate data
% lat_data =  wfield( mask, nsubj );
% 
% % Convolution fields
% cfield  = convfield( lat_data, params );
% dcfield = convfield( lat_data, params, 1 );
% 
% 
% % Regression estimators
% % values for regression
% uvals = -4:0.5:4;
% 
% % Which LKC needs to be estimated?
% L = [1, NaN, NaN];
% 
% L_reg_OLS  = LKC_regress_est( cfield, uvals, "OLS", L, 1, 1 );
% L_breg_OLS = LKC_regress_est( cfield, uvals, "OLS", L, 1, 200 );
% L_reg_SD   = LKC_regress_est( cfield, uvals, "SD", L, 1, 1 );
% L_breg_SD  = LKC_regress_est( cfield, uvals, "SD", L, 1, 200 );
% L_reg_SC   = LKC_regress_est( cfield, uvals, "SC", L, 1, 1 );
% L_breg_SC  = LKC_regress_est( cfield, uvals, "SC", L, 1, 200 );
% L_reg_PI   = LKC_regress_est( cfield, uvals, "PI", L, 1, 1 );
% L_breg_PI  = LKC_regress_est( cfield, uvals, "PI", L, 1, 200 );
% 
% struct(...
%     'theoryL', theoryLKC,...
%     'L_reg_OLS', L_reg_OLS,...
%     'L_breg_OLS', L_reg_OLS,...
%     'L_reg_SD', L_reg_SD,...
%     'L_breg_SD', L_breg_SD,...
%     'L_reg_SC', L_reg_SC,...
%     'L_breg_SC', L_breg_SC,...
%     'L_reg_PI', L_reg_PI,...
%     'L_breg_PI', L_breg_PI...
% )
% 
% 
% % Regression estimators
% % values for regression
% uvals = [-5, 0, 5];
% 
% % Which LKC needs to be estimated?
% L = [1, NaN, theoryLKC(2)];
% 
% L_reg_OLS  = LKC_regress_est( cfield, uvals, "OLS", L, 1, 1 );
% L_breg_OLS = LKC_regress_est( cfield, uvals, "OLS", L, 1, 200 );
% L_reg_SD   = LKC_regress_est( cfield, uvals, "SD", L, 1, 1 );
% L_breg_SD  = LKC_regress_est( cfield, uvals, "SD", L, 1, 200 );
% L_reg_SC   = LKC_regress_est( cfield, uvals, "SC", L, 1, 1 );
% L_breg_SC  = LKC_regress_est( cfield, uvals, "SC", L, 1, 200 );
% L_reg_PI   = LKC_regress_est( cfield, uvals, "PI", L, 1, 1 );
% L_breg_PI  = LKC_regress_est( cfield, uvals, "PI", L, 1, 200 );
% 
% struct(...
%     'theoryL', theoryLKC,...
%     'L_reg_OLS', L_reg_OLS,...
%     'L_breg_OLS', L_reg_OLS,...
%     'L_reg_SD', L_reg_SD,...
%     'L_breg_SD', L_breg_SD,...
%     'L_reg_SC', L_reg_SC,...
%     'L_breg_SC', L_breg_SC,...
%     'L_reg_PI', L_reg_PI,...
%     'L_breg_PI', L_breg_PI...
% )
% -------------------------------------------------------------------------
% AUTHOR: Fabian Telschow
%--------------------------------------------------------------------------
%% Add and check optional values
%--------------------------------------------------------------------------
% Get constants from the input
sY = field.fieldsize;
D  = field.D;
N  = field.fibersize;

if ~exist( 'uvals', 'var' )
   % default number of bootstrap replicates
   uvals = -1:0.2:1;
end

if ~exist( 'CovEst', 'var' )
   % default number of bootstrap replicates
   CovEst = "OLS";
end

if ~exist( 'L', 'var' )
   % default number of bootstrap replicates
   L = [1, NaN(1, D)];
end

if ~exist( 'Mboot', 'var' )
   % default number of bootstrap replicates
   Mboot = 1;
end

if N == 1 && Mboot > 1
    error('The bootstrap version requires that N > 2!')
end

if ~exist( 'standardize', 'var' )
    % Default value of "normalize"
    if N == 1
       standardize = 0;
    else
       standardize = 1;
    end
end

if ~exist( 'version', 'var' )
    % Default value of "version"
    version = "C";
end

%% Check mandatory input and get important constants
%--------------------------------------------------------------------------

lU = length(uvals);

% Check that method is implemented for dimension D
if D > 3
    error( 'D must be < 4. Higher dimensional domains have not been implemented')
end


%% Main function
%--------------------------------------------------------------------------
indexD = repmat(":", [1, D]);

% Mask the field
field = Mask(field, -Inf);

% Get the Empirical estimate of the expected Euler charactersitic either
% through a Gaussian bootstrap or from the fields directly
if Mboot > 1
    % Save the field as a matrix
    field = field.field;
    
    % Get weights for the multiplier bootstrap
    multiplier = normrnd( 0, 1, [ N, Mboot ] );
    
    % Reshape and and standardize the field, such that it has unit variance
    field = reshape( field, prod( sY( 1:end-1 ) ), N );
    
    % Normalize the residuals
    field = field - mean( field, 2 );
    field = field ./ sqrt( sum( field.^2, 2 ) );
    field(isnan(field)) = -Inf;

  %  ff = NaN([sY(1:2), Mboot]);
    
    for i = 1:Mboot
        % Get the bootstrapped process
        if D ~= 1
            mfield = reshape( field * multiplier( :, i ), sY( 1:end-1 ) );
        else
            mfield = field * multiplier( :, i );
        end
        
        if i == 1
            EC = ECcurve( Field(mfield, D), uvals);
        else
            ECi = ECcurve( Field(mfield, D), uvals);
            EC.field = [EC.field, ECi.field];
        end
        
         %   ff(:,:,i) = mfield;
    end
   
    % Indicate that the number of regression curves changed 
    N = Mboot;
    
else
    % Normalize the residuals
    if standardize
        field = normalize_field(field);
    end

    % make sure NaN's are -inf
    field.field(isnan(field.field)) = -Inf;
    
    for i = 1:N    
        if i == 1
            EC = ECcurve(field(indexD{:}, i), uvals);
        else
            ECi = ECcurve(field(indexD{:}, i), uvals);
            EC.field = [EC.field, ECi.field];
        end
    end
end

% Subtract the known LKCs from the EC curve
L_exist = find(~isnan(L));
L_empty = find(isnan(L));

if ~isempty(L_exist)
    for d = L_exist
        rho_d = EC_denstity( d-1 );
        
        EC.field = EC.field - L(d) * rho_d(uvals)';
    end
end

% Get the mean Euler characteristic for regression!
mEC = mean(EC);

% Get the design matrix
if ~isempty(L_empty)
    X = [];
    for d = L_empty
        rho_d = EC_denstity( d-1 );
        X     = [X, rho_d(uvals)'];
    end
end

% Get the estimate of the Covariance
switch CovEst
    case "OLS" % Ordinary least square
        C1 = eye(lU);
    case "SD" % smoothed diagonal
        params = ConvFieldParams( 1, 0 );
        sv = convfield( var(EC), params );
        C1 = diag(1 ./ sv.field  );
        C1(isinf(C1)) = 0;
    case "SC"
        % Smooth the empirical covariance using a Gaussian kernel
        params    = ConvFieldParams( [1 1], 0 );
        xvals     = cell([1 2]);
        xvals{1}  = EC.xvals{1};
        xvals{2}  = EC.xvals{1};
        C1 = convfield(  Field( cov(EC.field'), xvals ), params );
        C1 = pinv( C1.field );
        clear SC xvals
    case "PI"
        C1 = pinv(cov(EC.field') );
end

% Get the GLS estimate
XCX1    = pinv(X' * C1 * X);
LKC_est = XCX1 * X' * C1 * mEC.field;

% Get the estimated variance of the estimator
R = EC.field - X * LKC_est;
s2_est = mean(trace(R' * C1 * R)) / (lU - rank(X) - 1 ); % check degrees of freedom formula in GLS!
LKC_est_var = s2_est * XCX1;

L(L_empty) = LKC_est';
LKC_est = L(2:end);
return