function [ inv ] = invPerm( perm )
% invperm cacluates the inverse of a permutation.
%--------------------------------------------------------------------------
% ARGUMENTS
% perm  A permutation to be inverted.
%--------------------------------------------------------------------------
% OUTPUT
% inv   the inverse permutation.
%--------------------------------------------------------------------------
% EXAMPLES
% invPerm([4,3,1,2])
% 
% a = [4,6,8,2]
% b = a([1,3,4,2])
% b(invPerm([1,3,4,2]))
%
% a = [4,6,8,2];
% [sorted_a,sort_index] = sort(a);
% sorted_a(invPerm(sort_index))
%--------------------------------------------------------------------------

n = length(perm);
inv = 1:n;
inv(perm) = 1:n;

end

