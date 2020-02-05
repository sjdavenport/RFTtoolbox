function [ lmarrayindices, lmInd ] = lmindices_new( Y, top, mask, CC )
% lmindices_new( Y, top, CC, mask ) finds the top local maxima of an image,
% given a conectivity criterion CC, that lie within a specified mask.
%--------------------------------------------------------------------------
% ARGUMENTS
% Y      a 3 dimensional array of real values
% top    the top number of local maxima of which to report
% CC     the connectivity criterion
% mask   a 0/1 mask.
%--------------------------------------------------------------------------
% OUTPUT
% lmarrayindices    an npeaks by D array of the peak locations
% lmInd 
%--------------------------------------------------------------------------
% EXAMPLES
% a = zeros([91,109,91]);
% a(16,100,40) = 5;
% a(10,50,35) = 3;
% lmindices_new(a,2)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
if nargin < 2
    top = 1;
end
if nargin < 3
    mask = ones(size(Y));
end
if nargin < 4
    CC = 26;
end

Ydim = size(Y);

r = 0.1;
[a,b,c] = ndgrid(r*(1:size(Y,1)), r*(1:size(Y,2)), r*(1:size(Y,3)));
ramp = min(Y(:)) - (a+b+c);
Y(logical(1-mask)) = ramp(logical(1-mask));

indices = find(imregionalmax(Y, CC).*mask); 
%ensures that the boundary of the mask are not calculated as local maxima 

[~, sorted_lm_indices] = sort(Y(indices), 'descend');
top = min(top, length(sorted_lm_indices));
lmInd = indices(sorted_lm_indices(1:top)); %choose the top indices and return with top first second second etcetera

[lmarrayindicesx,lmarrayindicesy, lmarrayindicesz]  = ind2sub(Ydim, lmInd);
lmarrayindices = [lmarrayindicesx, lmarrayindicesy, lmarrayindicesz];

end

