function vec = zero2nan( vec )
% ZERO2NAN(vec) sets all of the zeros in vec to NaN
%--------------------------------------------------------------------------
% ARGUMENTS
% vec       An array of numbers
%--------------------------------------------------------------------------
% OUTPUT
% vec       vec converted so that all the zeros are set to NaN
%--------------------------------------------------------------------------
% EXAMPLES
% zero2nan([1,5.1,0, 0, 0.0001])
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

if islogical(vec)
    vec = vec + 0; % Convert logical array to numeric (is there a better way?)
end

vec(abs(vec)<=2*eps) = nan;

end

