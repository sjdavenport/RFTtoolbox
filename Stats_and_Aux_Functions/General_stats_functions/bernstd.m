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
% bernstd(0.05,1000,0.95)
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

std_error = zeros(1,length(p));
interval = zeros(2,length(p));

%%  Main Function Loop
%--------------------------------------------------------------------------
for I = 1:length(p)
    std_error(I) = (p(I)*(1-p(I)))^(1/2)*norminv( 1-(1-level)/2 )/sqrt(N);
    interval(:,I) = [ p(I) - std_error(I), p(I) + std_error(I) ];
end
end

