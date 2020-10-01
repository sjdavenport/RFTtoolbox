function [ interval, std_error ] = bernstd( p, N, alpha )
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
if ~exist( 'alpha', 'var' )
   % default option of opt1
   alpha = 0.95;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
std_error = (p*(1-p))^(1/2)*norminv( 1-(1-alpha)/2 )/sqrt(N);
interval = [ p - std_error, p + std_error ];

end

