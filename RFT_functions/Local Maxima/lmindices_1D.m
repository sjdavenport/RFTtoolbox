function [ lmInd, IntInd ] = lmindices_1D(Y, top, mask)
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
if size(Y, 2) == 1
    Y = Y';
end

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
[~, Yindices] = sort(Y(indices), 'descend'); 

top = min(top, length(local_maxima));
lmInd = sort(local_maxima(Yindices(1:top)));

potentialBInds = lmInd;
BInd = [];
if lmInd(1) == 1
    BInd(1) = lmInd(1);
    potentialBInds = potentialBInds(2:end);
end
if lmInd(end) == length(Y)
    BIndend = lmInd(end);
    potentialBInds = potentialBInds(1:end-1);
else
    BIndend = [];
end
for I = 1:length(potentialBInds)
    if isnan(Y(potentialBInds(I)-1)*Y(potentialBInds(I)+1))
        BInd = [BInd,potentialBInds(I)];
    end
end
BInd = [BInd, BIndend];
IntInd = setdiff(lmInd, BInd);

% [~, sorted_lm_indices] = sort(Y(indices), 'descend');
% top = min(top, length(sorted_lm_indices));
% lmInd = sorted_lm_indices(1:top);

end