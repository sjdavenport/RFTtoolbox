function [u_FWER, u_FWER_se] = get_EECthreshold( FWER, uvals, LKC, D, LKC0,...
                                                 type )
% get_EECthreshold( FWER, uvals, LKC, D, LKC0, type )
% calculates the EECthreshold.
%--------------------------------------------------------------------------
% ARGUMENTS
% FWER  an number between [0 1] controlling the FWER
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
% u_FWER  threshold controlling the FWER using the EEC as approximation for
%         the max distribution
% u_FWER_se  standard error of u_FWER threshold assuming asymptotic
%            normality.
%--------------------------------------------------------------------------
% EXAMPLES
% dim   = [50 50 50];
% nsubj = 30;
% D     = length( dim );
% uvals = -6 : 0.1 : 6;
% mask  = true(dim);
% 
% Y = noisegen( dim, nsubj, 5*2*sqrt(2*log(2)), 0 );
% HPE  = LKCestim_HPE( Y, D, mask, 1, "C" );
% 
% [u_FWER, u_FWER_se] = get_EECthreshold( 0.05, uvals, HPE.hat1, D, HPE.L0,...
%                                                  "gaussian" );
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow

%%%%%%%%%%%%%%%%%%%%%%%%%%%% add default values %%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5
    LKC0   = 1;
    type = "gaussian";
end
if nargin < 6
    type = "gaussian";
end

%%%%%%%%%%%%%% constants from input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sLKC = size(LKC);
N    = sLKC(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% main function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(type, "gaussian")
    [ ~, hatEEC, ~, C ] = EEC( uvals, LKC, D, LKC0, type );
    [ind, u0]           = crossing( hatEEC, uvals, FWER, 'linear' );
    u_FWER              = max(u0);
end

% compute the variability of the threshold
if N > 1
    tau       = [ 0; diff( hatEEC ) ] ./ diff(uvals);
    u_FWER_se = sqrt( C( max(ind), max(ind ) ) ...
                            / ( N * tau(max(ind))^2 ) );
else
    u_FWER_se = NaN;
end