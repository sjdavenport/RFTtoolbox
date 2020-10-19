function [max_coverage, cluster_coverage, x] = ECcurveanal( coverage, sample_size, L, L0, quantile, top )
% ECcurveanal( coverage, mask, sample_size, quantile, L ) plots the upper
% tail of the EC curves.
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
if ~exist( 'quantile', 'var' )
    % default option of opt1
    quantile = 0.1;
end

if ~exist( 'L', 'var')
    L = mean(coverage.storeLKCs,2)';
end

if ~exist( 'top', 'var')
    top = size(coverage.allmaxima,1);
end

%%  Main Function Loop
%--------------------------------------------------------------------------

if quantile > 0
    [ curve, x ] = maxECcurve( coverage.convmaxima, quantile );
    [ curve_all, x_all ] = maxECcurve( coverage.allmaxima, quantile );
    
    plot(x, curve);
    curve_conv = EEC( x, L, L0, 'T', sample_size -1 );
    hold on
    plot(x,curve_conv)
    hold on
    plot(x_all(end-length(x):end),curve_all(end-length(x):end))
    legend('max', 'LKC', 'all')
end

threshold = EECthreshold( 0.05, L, L0, 'T', sample_size -1)
threshold_p = EECthreshold( 0.05, L, L0, 'T', sample_size -1, 'p')

if top < size(coverage.allmaxima,1)
    coverage.allmaxima = coverage.allmaxima(1:top,:);
end
max_coverage = sum(coverage.convmaxima > threshold)/size(coverage.convmaxima, 2)
cluster_coverage = sum(coverage.allmaxima(:) > threshold)/(size(coverage.allmaxima, 2))

end

