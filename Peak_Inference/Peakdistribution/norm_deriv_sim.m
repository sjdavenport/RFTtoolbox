function [ field_sim_store, deriv2_sim_store ] = norm_deriv_sim( omega, Lambda, sigma, niters )
% NEWFUN serves as a function template.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'niters', 'var' )
   % Default value
   niters = 1000;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
field_sim_store = normrnd(0,sigma,1,niters);
deriv2_sim_store = -field_sim_store*Lambda +  normrnd(0, sqrt(omega - Lambda^2/sigma^2), 1, niters);

end

