function [ lmarrayindices, lmind, peakvals ] = lmindices( Y, top, mask, CC )
% lmindices( Y, top, CC, mask ) finds the top local maxima of an image,
% given a conectivity criterion CC, that lie within a specified mask.
% NEED TO REWRITE SO THAT IT TAKES IN A FIELD OBJECT!!!
%--------------------------------------------------------------------------
% ARGUMENTS
% Y      a 3 dimensional array of real values
% top    the top number of local maxima of which to report
% CC     the connectivity criterion
% mask   a 0/1 mask
%--------------------------------------------------------------------------
% OUTPUT
% lmarrayindices    an D by npeaks array of the peak locations
% lmind             an npeaks by 1 array where each value is the converted
%                   peak location, converted according to the size of Y 
%                   using convind
% peakvals          an npeaks by 1 array each value of which is the value
%                   the field takes at a given peak
%--------------------------------------------------------------------------
% EXAMPLES
% %3D example
% a = zeros([91,109,91]);
% a(16,100,40) = 5;
% a(10,50,35) = 3;
% [peaklocs, peakinds, peakvals] = lmindices(a,2)
% a(peakinds) %== peakvals
%
% %1D example 
% lmindices([1,2,1])
% lmindices([1,2,1,2,1],2)
%
% %2D example
% [ lmarrayindices, lmInd ] = lmindices([1,1,1;1,2,1;1,1,1])
% lmindices([1,1,1;1,2,1;1,1,1;1,1,2], 2)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if ~exist('top', 'var')
    top = 1;
elseif strcmp(top, 'all')
    top = numel(Y);
end
Y = squeeze(Y);
Ydim = size(Y);
D = length(Ydim); %Note in the D = 1 case we still take D = 2 and use a 
%connectivity criterion of 8 as this is equivalent to a connectivity criterion of 2!

if ~exist('mask', 'var')
    mask = ones(Ydim);
end
if ~isequal(size(mask), size(Y))
    error('The mask must be the same size as Y')
end

if ~exist('CC', 'var')
    if D == 2
        CC = 8;
    elseif D == 3
        CC = 26;
    else
        error('NotworkedoutD>3yet')
    end
end

r = 0.1;
[a,b,c] = ndgrid(r*(1:size(Y,1)), r*(1:size(Y,2)), r*(1:size(Y,3)));
ramp = min(Y(:)) - (a+b+c);
Y(logical(1-mask)) = ramp(logical(1-mask));

indices = find(imregionalmax(Y, CC).*mask); 
%ensures that the boundary of the mask are not calculated as local maxima 

[~, sorted_lm_indices] = sort(Y(indices), 'descend');
top = min(top, length(sorted_lm_indices));
lmind = indices(sorted_lm_indices(1:top)); %choose the top indices and return with top first second second etcetera

if D == 2
    [lmarrayindicesx,lmarrayindicesy]  = ind2sub(Ydim, lmind);
    lmarrayindices = [lmarrayindicesx, lmarrayindicesy]';
elseif D == 3
    [lmarrayindicesx,lmarrayindicesy, lmarrayindicesz]  = ind2sub(Ydim, lmind);
    lmarrayindices = [lmarrayindicesx, lmarrayindicesy, lmarrayindicesz]';
end

if D == 2
    if Ydim(1) == 1
        lmarrayindices = lmarrayindicesy;
    elseif Ydim(2) == 1
        lmarrayindices = lmarrayindicesx;
    end
end

% Obtain the values at the peaks
peakvals = Y(lmind);

% Return a column vector if the dimension is 1
% if D == 1
%     lmarrayindices = lmarrayindices';
%     lmind = lmind';
%     peakvals = peakvals';
% end

end
