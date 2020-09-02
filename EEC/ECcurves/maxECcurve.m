function [ curve, x ] = maxECcurve( maxima, alpha_percentile, increm )
% maxECcurve( maxima ) estimates the tail of the EC curve from the
% distribution of observed maxima over many iid iterations
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  maxima
%  alpha_percentile
% Optional
%  increm
%--------------------------------------------------------------------------
% OUTPUT
%  curve     the Euler characteristic curve
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

curve = zeros(1, length(x));
for I = 1:length(x)
    curve(I) = sum(maxima(:)>x(I));
end
curve = curve/n_maxima;

end

