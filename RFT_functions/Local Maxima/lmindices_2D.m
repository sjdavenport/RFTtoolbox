function [lmarrayindices, lmInd] = lmindices_2D(Y, top, CC, mask)
% lmindices_2D(Y, top, CC, mask)
%--------------------------------------------------------------------------
% ARGUMENTS
% Y      a 2 dimensional array of real values
% top    the top number of local maxima of which to report
% CC     the connectivity criterion
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
    CC = 4;
end
if nargin < 4
    mask = ones(size(Y));
end
mask = zero2nan(mask);

Y = Y.*mask;
Y(isnan(Y)) = min(Y(~isnan(Y))) - 1;

Ydim = size(Y);

dYx = diff(Y, 1, 2); %difference is taken in the x direction i.e. horizontal
dYy = diff(Y, 1, 1); %difference is taken in the y direction i.e. vertical
dYxpos = (dYx > 0);
dYxneg = (dYx < 0);
dYypos = (dYy > 0);
dYyneg = (dYy < 0);

rowwiselms = [ones(Ydim(1), 1), dYxpos].*[dYxneg, ones(Ydim(1), 1)];
colwiselms = [ones(1, Ydim(2)); dYypos].*[dYyneg; ones(1, Ydim(2))];

if CC == 8
    error('Not set up yet')
    %something like this is need:
    %     dYdiagse = Y[2:Ydim[1], 2:Ydim[2]] - Y[1:(Ydim[1]-1), 1:(Ydim[2]-1)]
    %     dYdiagnw = -dYdiagse
    %
    %     diagwise = rbind(matrix(rep(TRUE, Ydim[2]), 1, Ydim[2]), cbind(matrix(rep(TRUE, Ydim[1] - 1), Ydim[1] - 1), dYdiagse))
    %     diagwise = diagwise*rbind(cbind(dYdiagnw, matrix(rep(TRUE, Ydim[1] - 1), Ydim[1] - 1)), matrix(rep(TRUE, Ydim[2]), 1, Ydim[2]))
    %
    %     dYdiagne = Y[1:(Ydim[1]-1), 2:Ydim[2]] - Y[2:Ydim[1], 1:(Ydim[2]-1)]
    %     dYdiagsw = -dYdiagne
    %
    %     diagwise = diagwise*rbind(cbind(matrix(rep(TRUE, Ydim[1] - 1), Ydim[1] - 1), dYdiagne), matrix(rep(TRUE, Ydim[2]), 1, Ydim[2]))
    %     diagwise = diagwise*rbind(matrix(rep(TRUE, Ydim[2]), 1, Ydim[2]), cbind(dYdiagsw, matrix(rep(TRUE, Ydim[1] - 1), Ydim[1] - 1)))
    %
    %     indices = which( rowwiselms*colwiselms*diagwise > 0 )
    %     array_indices = which( rowwiselms*colwiselms*diagwise > 0, arr.ind = T )
else
    indices = find( rowwiselms.*colwiselms > 0 );
end

% [array_indices_row, array_indices_col]  = ind2sub(Ydim, indices);

[~, sorted_lm_indices] = sort(Y(indices), 'descend');
top = min(top, length(sorted_lm_indices));
lmInd = indices(sorted_lm_indices(1:top)); %choose the top indices and return with top first second second etcetera

[lmarrayindicesrow,lmarrayindicescol]  = ind2sub(Ydim, lmInd);
lmarrayindices = [lmarrayindicesrow, lmarrayindicescol];
end