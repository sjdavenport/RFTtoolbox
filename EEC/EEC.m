function [ EEC_vals, mEEC, EEC_se_hat, C_hat ] = EEC( uvals, LKC, LKC0,...
                                                      type, df, outfield )
% EEC( uvals, LKC, D, LKC0, type )
% calculates the expected Euler characteristic curves for given LKC values
% at all values uvals
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   uvals	an Nuvals length vector whose entries are the real numbers, where
%           the EEC is evaluated
%   LKC     an D by N length vector containing N estimates of the Lipschitz
%           Killing curvatures.
%   LKC0    Euler characteristic of the domain
%
% Optional
%   type	string specifying the type of random field the EEC is computed from.
%           Default option is "gaussian". (later "t","F","") will follow.
%   df      value depends on type
%           - "t" an numeric the degrees of freedom.
%--------------------------------------------------------------------------
% OUTPUT
% EEC_vals    an Nuvals x N array containing the EEC curves evaluated
%             at the locations specified by uvals
% mEEC        mean of the EEC curves evaluated at uvals
% EEC_se_hat  standard error of mEEC evaluated at uvals
% C_hat       estimate of covariance function of the EC process
%--------------------------------------------------------------------------
% EXAMPLES
% dim   = [50 50 50];
% nsubj = 70;
% D     = length( dim );
% uvals = -6 : 0.1 : 6;
% mask  = true(dim);
% 
% Y    = noisegen( dim, nsubj, 5*2*sqrt(2*log(2)), 0 );
% HPE  = LKCestim_HPE( Y, mask, 1 );
% bHPE = LKCestim_HPE( Y, mask, 5e3 );
% 
% [ EEC_vals, mEEC, EEC_se_hat, C_hat ] = EEC( uvals, HPE.hatL1, D,...
%                                             HPE.L0, "gaussian" );
% 
% figure(1)                                         
% plot(uvals, EEC_vals, 'Color', ones([1 3])*0.75)
% line(uvals, mEEC', 'Color', 'blue', 'LineWidth', 1.5)
% line(uvals, mEEC'+ 1.96* [EEC_se_hat, -EEC_se_hat]', 'Color', 'blue',...
% 'LineWidth', 1.5, 'LineStyle', '--')
% title("HPE estimates of expected EC curve with pointwise CI's")
%
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow
%--------------------------------------------------------------------------

%% Check mandatory input and get important constants
%--------------------------------------------------------------------------

if ~exist( 'type', 'var' )
    type = "gaussian";
end

sLKC = size(LKC);
N    = sLKC(1);
D    = sLKC(2);

%% Check optional input
%--------------------------------------------------------------------------

if ~exist( 'outfield', 'var' )
    outfield = 1;
end

%% Main function
%--------------------------------------------------------------------------

% Create output field
EEC_vals = Field( { uvals } );

% Compute the EC curves
rho = [ ECdensity( uvals, 0, type, df );
        ECdensity( uvals, 1, type, df );...
        ECdensity( uvals, 2, type, df );...
        ECdensity( uvals, 3, type, df ) ]';

EEC_vals.field = rho(:, 1:(D+1)) * [ LKC0 * ones(N,1), LKC]';

% compute estimated se of EEC estimator and standard error
if N > 1
    % covariance of LKC estimator
    Sigma_hat = cov( LKC );
    % mean of EEC curve
    mEEC = mean( EEC_vals, 1 );
    % covariance matrix of EEC curve
    C_hat = rho( :, 1:D ) * Sigma_hat * rho( :, 1:D )';
    % standard error of EEC curve
    EEC_se_hat = sqrt( diag( C_hat ) / N );
else
    C_hat      = NaN;
    EEC_se_hat = NaN;
    mEEC       = EEC_vals;
end

if ~outfield
    EEC_vals = EEC_vals.field;
end

return