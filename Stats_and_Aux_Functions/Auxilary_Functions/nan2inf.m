function vec = nan2inf( vec, posorneg )
% INF2NAN(vec) sets all of the nans to zeros.
%--------------------------------------------------------------------------
% ARGUMENTS
% vec       A vector
%--------------------------------------------------------------------------
% OUTPUT
% vec       a vector whose infinities have been converted to NaNs
%--------------------------------------------------------------------------
% EXAMPLES
% nan2inf([1,5.1,5, NaN, 0.0001], -1)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if ~exist('posorneg', 'var')
    posorneg = 1;
end

vec(isnan(vec)) = posorneg*Inf;

end

