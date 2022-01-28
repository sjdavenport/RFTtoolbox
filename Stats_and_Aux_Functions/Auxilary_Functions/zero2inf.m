function vec = zero2inf( vec, posorneg )
% INF2NAN(vec) sets all of the nans to zeros.
%--------------------------------------------------------------------------
% ARGUMENTS
% vec       A vector
% posorneg  Specifies whether to use positive or negative infinity
%--------------------------------------------------------------------------
% OUTPUT
% vec       a vector whose infinities have been converted to NaNs
%--------------------------------------------------------------------------
% EXAMPLES
% zero2inf([1,5.1,5, 0, 0.0001], -1)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if ~exist('posorneg', 'var')
    posorneg = 1;
end

vec(vec == 0) = posorneg*Inf;

end

