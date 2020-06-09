function vec = inf2nan( vec )
% INF2NAN(vec) sets all of the nans to zeros.
%--------------------------------------------------------------------------
% ARGUMENTS
% vec       A vector
%--------------------------------------------------------------------------
% OUTPUT
% vec       a vector whose infinities have been converted to NaNs
%--------------------------------------------------------------------------
% EXAMPLES
% inf2nan([1,5.1,5, Inf, 0.0001])
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
vec(isinf(vec)) = NaN;

end

