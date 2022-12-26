function out = vec2lim( vec )
% VEC2LIM( vec ) takes a vector and returns a vector with the first and 
% last entries
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   vec      a vector of any length
%--------------------------------------------------------------------------
% OUTPUT
%   out      a vector with length 2, containing the first and last entries of vec
%--------------------------------------------------------------------------
% EXAMPLES
% vec2lim([1, 2, 3, 4]) % returns [1, 4]
% vec2lim([1, 1, 1, 1, 1]) % returns [1, 1]
% vec2lim([]) % returns []
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

out = [vec(1), vec(end)];

end

