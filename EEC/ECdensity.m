function rho_d = ECdensity( uvals, d, type, df )
% NEWFUN( var1, var2, opt1 ) is described here [...]
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   var1   this is a mandatory variable.
%   var2   this is a mandatory variable.
%
% Optional
%   opt1   this is an optional parameter. Default 0.
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

% This kind of code with exists is better than using nargin < xy, since
% Parameter can be easily permuted
if ~exist( 'type', 'var' )
   % Default option of opt1
   type = "gaussian";
   df   = 0;
end


%% Main function  
%--------------------------------------------------------------------------

% Compute the EC curves
if strcmp( type, "gaussian" )
    switch d
        case 0
            rho_d = 1 - normcdf( uvals );
        case 1
            rho_d = ones( 1, length( uvals ) ) .* ...
                    exp( -uvals.^2/2 ) / (2*pi)^(2/2);
        case 2
            rho_d = uvals .* exp( -uvals.^2/2 ) / (2*pi)^(3/2);
        case 3
            rho_d = ( uvals.^2 - 1 ) .* exp( -uvals.^2/2 ) / (2*pi)^(4/2);
    end
    
elseif strcmp( type, "t" )
    % Factor depending on degrees of freedom
    if df > 250
        fac = 1;
    else
        fac = gamma( (df + 1) / 2 ) / gamma( df / 2 ) / sqrt(df/2);
    end
    
    % Compute the ECdensities
    switch d
        case 0
            rho_d = tcdf( uvals, df, 'upper' );
        case 1
            rho_d = ones( 1, length( uvals ) ) .* ...
                    ( 1 + uvals.^2 / df ).^( - ( df - 1 ) / 2 ) / (2*pi)^(2/2);
        case 2
            rho_d = fac * uvals .*...
                    ( 1 + uvals.^2 / df ).^( - ( df - 1 ) / 2 ) / (2*pi)^(3/2);
        case 3
            rho_d = ( (df-1) / df * uvals.^2 - 1 ) .*...
                    ( 1 + uvals.^2 / df ).^( - ( df - 1 ) / 2 ) / (2*pi)^(4/2);
    end
end

return