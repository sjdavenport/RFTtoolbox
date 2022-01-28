function KRexpectation = KRdensityall( fV, dV, d2V, fd2cov, fdcov, dd2cov, u, niters )
% KRDENSITY() evaluates the density from the Kac-Rice formula.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Lambda    the covariance of the derivative
% Delta     the covariance of the second derivative
% gamma     the covariance between the field and the second derivative
% omega     the covariance between the first and second derivatives
% u         a vector giving the thresholds 
% niters    the number of iterations used to calculate the expectation
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

% Set the correlation between the first and second derivative to be zero 
% (this is e.g. true under stationarity).
if ~exist('fdcov', 'var')
    fdcov = 0;
end

if ~exist('fdcov', 'var')
    dd2cov = 0;
end

if ~exist('niters', 'var')
    niters = 1000000;
end

% If u is not specified just compute the epected number of local maxima
% without any restriction on their height, this is the same as taking u=-Inf
if ~exist('u', 'var')
    u = -Inf;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
covmate = [fV,fd2cov;fd2cov, d2V] - [fdcov;dd2cov]*[fdcov,dd2cov]/dV;
data = mvnrnd([0;0], covmate, niters);

setbelowzero = data(:,2) < 0;
KRexpectation = zeros(1, length(u));

for I = 1:length(u)
    setaboveu = data(:,1) > u(I);
    set2use = setbelowzero.*setaboveu;

    if sum(set2use) == 0
        KRexpectation(I) = 0;
    else
        temp_data = data(:,2);
        temp_data(~set2use) = 0;
        KRexpectation(I) = mean(abs(temp_data));
    end
end

% if sum(setbelowzero) == 0
%     KRexpectation = 0;
% else
%     KRexpectation = mean(abs(data(setbelowzero)));
% end

KRexpectation = KRexpectation/sqrt(2*pi*dV); %multiply by p_t(0)!

end
