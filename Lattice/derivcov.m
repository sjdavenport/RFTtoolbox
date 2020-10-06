function [ Delta, Lambda, Omega ] = derivcov( fields, lag )
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
if ~exist( 'opt1', 'var' )
   % default option of opt1
   opt1 = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
N = size(fields,2);
deriv = diff(fields,1);
deriv2 = diff(deriv,1);

deriv = deriv - mean(deriv,2);
Lambda = var(deriv,0,2);
Omega = var(deriv2,0,2);

if lag > 0
    Delta = mean(deriv(1:end-lag,:).*deriv2(lag:end, :), 2)*(N-1)/N;
else
    lag = -lag;
    Delta = mean(deriv((lag+1):end,:).*deriv2(1:end-lag+1,:), 2)*(N-1)/N;
end

end

