function u_FWER = EECthreshold( FWER, LKC, LKC0, type, df, option )
% EECthreshold( FWER, LKC, LKC0, type, df ) this function computes the FWER
% threshold from the EC heuristic.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   FWER  an numeric strictly between 0 and 1 denoting the FWER which
%         needs to be controlled.
%   LKC   an 1 x D vector containing the first to dth LKCs.
%   LKC0  an integer containing the zeroth LKC
%
% Optional
%   type   a character indicating the type of EEC curve used for the
%          threshold. Options are same as in EEC(). Default "Z".
%   df     degrees of freedom for the chosen type, see EEC().
%
%--------------------------------------------------------------------------
% OUTPUT
%   u_FWER  the FWER threshold
%
%--------------------------------------------------------------------------
% DEVELOPER TODOs:
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow  
%--------------------------------------------------------------------------


%% Check mandatory input and get important constants
%--------------------------------------------------------------------------


%% Add/check optional values
%--------------------------------------------------------------------------

if ~exist( 'type', 'var' )
   % Default option of type
   type = "Z";
end

if ~exist( 'df', 'var' )
   % Default option of type
   df = 2;
end

if ~exist( 'option', 'var' )
   option = 'standard'; 
end

%% Main function  
%--------------------------------------------------------------------------

if strcmp(option, 'standard')
    EECminusalpha = @(x) EEC( x, LKC, LKC0, type, df, 0 ) - FWER;
    u_FWER = fzero( EECminusalpha, [ 1, 200 ] );
elseif strcmp(option, 'poisson') || strcmp(option, 'Poisson') || strcmp(option, 'P') || strcmp(option, 'p')
    possionprob = @(x) 1-exp(-EEC( x, LKC, LKC0, type, df, 0 )) - FWER;
    u_FWER = fzero( possionprob, [ 2, 200 ] );
else 
    error('This option has not been set')
end
return