function rho_d = EC_denstity( d, type )
% rho_d = EC_denstity( d, type ) returns a function handle computing the
% EC density.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  d     d-th EC density. (Currently d<4)
%  type  type of the field. Currently only "Z"-field supported.
%
%--------------------------------------------------------------------------
% OUTPUT
%   rho_d a function handle computing the d-th EC density of the field
% -------------------------------------------------------------------------
% DEVELOPER TODOs:
% - write meaningful examples
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow
%--------------------------------------------------------------------------

%% Add and check optional values
%--------------------------------------------------------------------------

if ~exist( 'type', 'var' )
    % default option of limits
    type = "Z";
end

%% Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Check that d is smaller than 4
if d > 3
    error( 'd must be < 4. Higher dimensional domains have not been implemented.')
end

% Check that type is a string
if ~isstring(type)
    error('type must be a string indicating the type of the field. Possible values are Z, ...')
end

%% Main function
%--------------------------------------------------------------------------

if type == "Z" 
    % get the EC density
    switch d
        case 0
            rho_d = @(x) 1 - normcdf(x);
        case 1
            rho_d = @(x) (2*pi)^( -( d + 1 ) / 2 ) .* exp( -x.^2 / 2 );
        case 2
            rho_d = @(x) x * (2*pi)^( -( d + 1 ) / 2 ) .* exp( -x.^2 / 2 );
        case 3
            rho_d = @(x) (x.^2 - 1) * (2*pi)^( -( d + 1 ) / 2 ) .* exp( -x.^2 / 2 );
    end
end

return