function [ meanECcurve, stdECcurve, x, set_of_curves] = tailECcurve( maxima, alpha_percentile, increm )
% TAILECCURVE( maxima, alpha_percentile, increm ) uses the maxima recorded 
% over many simulations to generate an empirical estimate of the tail of
% the Euler characteristic curve
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  maxima   the set of recorded maxima with which to construct the tail of 
%           the EC curve
% Optional 
%  alpha_percentile     the proportion at which to start from
%  increm   the increment at which to calculate the curves
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
if ~exist( 'increm', 'var' )
   % default option of increm
   increm = 0.01;
end

if ~exist( 'alpha_percentile', 'var')
    alpha_percentile = 0.05;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
limits = [prctile(maxima(:), 100*(1-alpha_percentile) ), max(maxima(:))];

x = limits(1):increm:limits(2);
n_maxima = size(maxima,2);

set_of_curves = zeros(length(x),n_maxima);
for I = 1:length(x)
    set_of_curves(I, :) = sum(maxima > x(I), 1);
end

meanECcurve = mean(set_of_curves, 2);
stdECcurve = std(set_of_curves, 0, 2);

end

