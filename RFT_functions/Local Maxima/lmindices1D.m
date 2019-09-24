function [ indices, npeaks ] = lmindices( map, top, mask )
% LMINDICES( map, top ) finds the indices of the top local maxima of
% a given map. It uses a connectivity criterion of 18.
%--------------------------------------------------------------------------
% ARGUMENTS
% map   a 3D map which could be vectorized (ie map = map(:)) but not
%       necessary I guess.
% top  specifies the number of the top local maxima to find if its
%       intially set as inf it then counts the number of peaks within the mask!
% mask  an image of 1s and 0s that masks the data. (Needs to be 91 by 109
%       by 91).
%--------------------------------------------------------------------------
% OUTPUT
% indices   a vector of indices of local maxima.
% npeaks    the number of peaks returned = length(indices) this is only
%           really relevant when top = Inf as then the number of peaks that
%           is returned is variable.
%--------------------------------------------------------------------------
% EXAMPLES

dY = diff(map)
dYpos = (dY > 0)
dYneg = (dY < 0)

indices = which( dYpos[1:length(dYpos)-1]*dYneg[2:length(dYneg)] > 0 ) + 1
Yindices = sort(Y[indices], decreasing = T, index.return = T)


indices = sort(indices[Yindices$ix[1:top]], decreasing = F)

end
