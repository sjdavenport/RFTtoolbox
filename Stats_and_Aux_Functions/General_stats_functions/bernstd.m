function [ interval, std_error ] = bernstd( p, N, level )
% BERNSTD( p, N, sigma ) generates the bernouilli standard error confidence
% intervals that arise from the central limit theorem.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% bernstd(0.95,50,1)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'level', 'var' )
   % default option of alpha
   level = 0.95;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
std_error = (p*(1-p))^(1/2)*norminv( 1-(1-level)/2 )/sqrt(N);
interval = [ p - std_error, p + std_error ];

end

