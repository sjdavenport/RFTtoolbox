function out = sbasis( d, D )
% SBASIS( d, D ) returns the dth D-dimensional standard basis vector
%--------------------------------------------------------------------------
% ARGUMENTS
% d   the basis vector to return
% D   the dimension
%--------------------------------------------------------------------------
% OUTPUT
% out   the standard basis vector
%--------------------------------------------------------------------------
% EXAMPLES
% sbasis(2, 4)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
out = zeros(D,1);
out(d) = 1;
end

