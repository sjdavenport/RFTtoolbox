function max_dist = simthresh( df, FWHM, niters, STAT )
% SIMTHRESH(groupsize, STAT, niters) estimates the height needed such that
% thresholdx100 percent of the time the map will lie above it.
%--------------------------------------------------------------------------
% ARGUMENTS
% df    If df == 1, a GRF is generated. If df = [1, n] a t field with
%       n degrees of freedom is generated. If df = [m, n] an F field 
%       with m and n degrees of freedom is generated.
% FWHM  The smoothness.
% niters    the number of random fields to generate.
% STAT      either 'Z', 'T', or 'F'
%--------------------------------------------------------------------------
% OUTPUT
% Height of the threshold.
%--------------------------------------------------------------------------
% EXAMPLES
% max_dist = simthresh( [1,4], 0, 10 )
%--------------------------------------------------------------------------
% AUTHOR: Sam Davenport.
if nargin < 4
    STAT = 'ZorT';
end

max_dist = zeros(1, niters);

for iter = 1:niters
    rf = genRF( 1, df, FWHM ); %This generates a t-field (Gaussian if df = 1).
    if strcmp(STAT, 'F')
        rf = rf.^2; %This squares the t-field to make it into an F-field.
        max_dist(iter) = max(rf(:));
    else
        max_dist(iter) = max(rf(:));
    end 
end

end

