function [lmarrayindices, lmInd] = lmindices_3D(Y, top, CC, mask)
% lmindices_3D(Y, top, CC, mask)
%--------------------------------------------------------------------------
% ARGUMENTS
% Y      a 3 dimensional array of real values
% top    the top number of local maxima of which to report
% CC     the connectivity criterion
% mask   a 0/1 mask.
%--------------------------------------------------------------------------
% OUTPUT
% array_indices
% indices
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Samuel J. Davenport
if nargin < 2
    top = 1;
end
if nargin < 3
    CC = 6;
end
if nargin < 4
    mask = ones(size(Y));
end
mask = zero2nan(mask);

Y = Y.*mask;
Y(isnan(Y)) = min(Y(~isnan(Y))) - 1;

Ydim = size(Y);

dYx = diff(Y, 1, 1); %difference is taken in the x direction
dYy = diff(Y, 1, 2); %difference is taken in the y direction
dYz = diff(Y, 1, 3); %difference is taken in the z direction
dYxpos = (dYx > 0);
dYxneg = (dYx < 0);
dYypos = (dYy > 0);
dYyneg = (dYy < 0);
dYzpos = (dYz > 0);
dYzneg = (dYz < 0);

xwiselms = cat(1, ones(1, Ydim(2), Ydim(3)), dYxpos).* cat(1, dYxneg, ones(1, Ydim(2), Ydim(3)));
ywiselms = cat(2, ones(Ydim(1), 1, Ydim(3)), dYypos).* cat(2, dYyneg, ones(Ydim(1), 1, Ydim(3)));
zwiselms = cat(3, ones(Ydim(1), Ydim(2), 1), dYzpos).* cat(3, dYzneg, ones(Ydim(1), Ydim(2), 1));

if CC ~= 6
    error('Not set for other CCs yet');
else
    indices = find( xwiselms.*ywiselms.*zwiselms > 0 );
end

% [array_indices_row, array_indices_col]  = ind2sub(Ydim, indices);

[~, sorted_lm_indices] = sort(Y(indices), 'descend');
top = min(top, length(sorted_lm_indices));
lmInd = indices(sorted_lm_indices(1:top)); %choose the top indices and return with top first second second etcetera

[lmarrayindicesx,lmarrayindicesy, lmarrayindicesz]  = ind2sub(Ydim, lmInd);
lmarrayindices = [lmarrayindicesx, lmarrayindicesy, lmarrayindicesz];
end