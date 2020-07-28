function u_FWER = EECthreshold( FWER, LKC, LKC0, type, df )
% EECthreshold( FWER, LKC, D, LKC0, type )
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   FWER  this is a mandatory variable.
%   LKC   this is a mandatory variable.
%   LKC0
%
% Optional
%   type   this is an optional parameter. Default "gaussian".
%   df     degrees of freedom if "t" or "F" type is chosen.
%
%--------------------------------------------------------------------------
% OUTPUT
%   out1  this is the first output
%   out2  this is the second output
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
   type = "gaussian";
end

if ~exist( 'df', 'var' )
   % Default option of type
   df = -100;
end

%% Main function  
%--------------------------------------------------------------------------

EECminusalpha = @(x) EEC( x, LKC, LKC0, type, df, 0 ) - FWER;

u_FWER = fzero( EECminusalpha, [ 2, 200 ] );

return