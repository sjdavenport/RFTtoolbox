function KRexpectation = KRdensity( muprime, muprime2, Lambda, Omega )
% KRDENSITY() evaluates the density from the Kac-Rice formula.
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


%%  Main Function Loop
%--------------------------------------------------------------------------
niters = 1000000;
data = randn([1,niters]);
data = data*sqrt(Omega) + muprime2;

setbelowzero = data < 0;
if sum(setbelowzero) == 0 
    KRexpectation = 0;
else
    data(~setbelowzero) = 0;
    KRexpectation = mean(abs(data));
end

% if sum(setbelowzero) == 0
%     KRexpectation = 0;
% else
%     KRexpectation = mean(abs(data(setbelowzero)));
% end

KRexpectation = KRexpectation*exp(-muprime^2/(2*Lambda))/sqrt(2*pi*Lambda); %multiply by p_t(0)!

end

