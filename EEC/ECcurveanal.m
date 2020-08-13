function L = ECcurveanal( coverage, mask, sample_size, quantile, L )
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
D = length(size(mask));
if D == 2
    if (size(mask,1) == 1) || (size(mask,2) == 1)
        D = 1;
    end
end


%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'quantile', 'var' )
   % default option of opt1
   quantile = 0.1;
end

if ~exist( 'L', 'var')
    L = mean(coverage.storeLKCs,2)';
end

%%  Main Function Loop
%--------------------------------------------------------------------------
L0 = EulerChar(mask, 0.5, D);
[ curve, x ] = maxECcurve( coverage.convmaxima, quantile );
[ curve_all, x_all ] = maxECcurve( coverage.allmaxima, quantile );

plot(x, curve);
curve_conv = EEC( x, L, L0, 'T', sample_size -1 );
hold on 
plot(x,curve_conv)
hold on 
plot(x_all(end-length(x):end),curve_all(end-length(x):end))
legend('max', 'LKC', 'all')

threshold = EECthreshold( 0.05, L, L0, 'T', sample_size -1)
threshold_p = EECthreshold( 0.05, L, L0, 'T', sample_size -1, 'p')

max_coverage = sum(coverage.convmaxima > threshold)/size(coverage.convmaxima, 2)
cluster_coverage = sum(coverage.allmaxima(:) > threshold)/(size(coverage.allmaxima, 2))

end

