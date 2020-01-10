 function [lmarrayindices, lmInd] = lmindices(Y, top, mask, CC)
%   find the local maxima in an N-dimensional
% array
%--------------------------------------------------------------------------
% Y      a D dimensional array of real values
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
% max_index = lmindices(a)
% 
% a(10,50,35) = 3;
% top_2_lms = lmindices(a, 2)
%
% Y = [1,1,1;1,2,1;1,1,1;1,2,1;1,1,1];
% lmindices_2D(Y,2)
%--------------------------------------------------------------------------
% AUTHOR: Samuel J. Davenport
dimY = size(Y);
nD = length(dimY);
if nargin < 2
    top = 1;
end
warning('if there is an issue here may have mask and CC the wrong way around')
if nargin < 3
    mask = ones(size(Y));
end

if nargin < 4
    if nD == 2
        CC = 4;
    elseif nD == 3
        CC = 6;
    end
end

if nD == 2
    if dimY(1) == 1 
        nD = 1;
    elseif dimY(2) == 1
        nD = 1;
        Y = Y';
    else 
        nD = 2;
    end
end

if nD == 1
    lmInd = lmindices_1D(Y, top, mask);
    lmarrayindices = lmInd; %arrayindices == indices in 1D!
elseif nD == 2
    [lmarrayindices, lmInd] = lmindices_2D(Y, top, CC, mask);
elseif nD == 3
    [lmarrayindices, lmInd] = lmindices_3D(Y, top, CC, mask);
else
    error('N > 3 is not supported yet')
end

% Need to add code to remove the edges here!
% if edges == 0
%     for I = 1:size(lmarrayindices, 1)
%         if 
%     end
% end

end

