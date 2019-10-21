function lmInd = lmindices_1D(Y, top, mask)
% lmindices_1D(Y, top, mask) find the local maxima in a 1D random lattice 
% field.
%--------------------------------------------------------------------------
% ARGUMENTS
% Y      a 1 dimensional array of real values
% top    the top number of local maxima of which to report
%--------------------------------------------------------------------------
% OUTPUT
% local_maxima      the indicies of the top local maxima of Y
%--------------------------------------------------------------------------
% EXAMPLES
% Y = [1,2,1, 2, 4, 3]
% a = lmindices_1D(Y)
% a = lmindices_1D(Y,2,[1,1,1,0,0,0])
% a = lmindices_1D(Y,2)
%--------------------------------------------------------------------------
% AUTHORS: Sam Davenport 

if nargin < 2
    top = 1;
end
if nargin < 3
    mask = ones(1,length(Y));
end
mask = zero2nan(mask);

Y = Y.*mask;
Y(isnan(Y)) = min(Y(~isnan(Y))) - 1;
    
dY = diff(Y);
dYpos = (dY > 0);
dYneg = (dY < 0);

indices = [1, dYpos].*[dYneg, 1] > 0;
local_maxima = find(indices);
% [~, Yindices] = sort(Y(indices), 'descend');

[~, sorted_lm_indices] = sort(Y(indices), 'descend');
top = min(top, length(sorted_lm_indices));
lmInd = indices(sorted_lm_indices(1:top));
% top = min(top, length(local_maxima));
% local_maxima = sort(local_maxima(Yindices(1:top)));

end