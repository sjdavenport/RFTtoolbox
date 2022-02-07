function [ U, deriv2_sim_store ] = chi2_deriv_sim( df, Lambda, MofLambda, niters )
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
U = chi2rnd(df-1,1, niters); %% changed from U = chi2rnd(df,1, niters)
H = normrnd(0, sqrt(MofLambda), 1, niters);
deriv2_sim_store = zeros(1, niters);
for I = 1:niters
%     modul(I, 100)
    %deriv_sim_store(I) = normrnd(0, 2*sqrt(U(I))*sqrt(Lambda));
    deriv2_sim_store(I) = 2*Lambda*chi2rnd(df-1) + 2*(- U(I)*Lambda) + 2*U(I)^(1/2)*H(I); %+ deriv_sim_store(I)^2/U(I)/2 
end

% KRexpectation = -mean(deriv2_sim_store(deriv2_sim_store < 0));

end

