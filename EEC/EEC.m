function [ EEC_vals, mEEC, EEC_se_hat, C_hat ] = EEC( uvals, LKC, D, LKC0, type )
% EEC( uvals, LKC, D, LKC0, type )
% calculates the expected Euler characteristic curves for given LKC values
% at all values uvals
%--------------------------------------------------------------------------
% ARGUMENTS
% uvals	an Nuvals length vector whose entries are the real numbers, where
%       the EEC is evaluated
% LKC   an D by N length vector containing N estimates of the Lipschitz
%       Killing curvatures.
% D     dimension of the domain
% LKC0  Euler characteristic of the domain
% type	string specifying the type of random field the EEC is computed from.
%       Default option is "gaussian". (later "t","F","") will follow.
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
%
%%%%%%%%%%%%%% add default values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    LKC0   = 1;
    type = "gaussian";
end
if nargin < 4
    type = "gaussian";
end

%%%%%%%%%%%%%% constants from input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sLKC = size(LKC);
N    = sLKC(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp( type, "gaussian" )
    % matrix containing the Hermite polynomials evaluated on uvals
    H      = [ ones( 1, length(uvals) ); uvals; ( uvals.^2 - 1 ) ]';
    % matrix containing the  EC densities evaluated on uvals
    rho    = H .* [ exp( -uvals.^2/2 ) / (2*pi)^(2/2);...
                    exp( -uvals.^2/2 ) / (2*pi)^(3/2);...
                    exp( -uvals.^2/2 ) / (2*pi)^(4/2) ]';
    % marginal part of GKF
    EEC_vals = ( 1 - normcdf( uvals ) )' * LKC0;
end

EEC_vals = EEC_vals + rho(:, 1:D) * LKC;

% compute estimated se of EEC estimator and standard error
if N > 1
    % covariance of LKC estimator
    Sigma_hat = cov( LKC' );
    % mean of EEC curve
    mEEC = mean( EEC_vals, 2 );
    % covariance matrix of EEC curve
    C_hat = rho( :, 1:D ) * Sigma_hat * rho( :, 1:D )';
    % standard error of EEC curve
    EEC_se_hat = sqrt( diag( C_hat ) / N );
else
    C_hat      = NaN;
    EEC_se_hat = NaN;
    mEEC       = EEC_vals;
end
return