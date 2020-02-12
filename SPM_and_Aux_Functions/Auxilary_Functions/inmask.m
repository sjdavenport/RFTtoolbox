function out = inmask(point, mask)
% INMASK(point, mask) tests whether a given point is in a given mask.
%--------------------------------------------------------------------------
% ARGUMENTS
% point
% mask
%--------------------------------------------------------------------------
% OUTPUT
% out
%--------------------------------------------------------------------------
% EXAMPLES
% %1D
% inmask( 0.9, ones(1,5))
% inmask( 0.45, ones(1,5))
% inmask( 5.45, ones(1,5))
% inmask( 5.55, [ones(1,5),0])
% 
% %2D
% mask = [0,0,1;0,0,1;0,0,1];
% inmask([1,2], mask)
% inmask([2.5,2.6], mask)
% inmask([2.5,2.4], mask)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
Dim = size(mask);
if Dim(1) == 1
    point = [1, point];
end
nearest_voxel = round(point);
if sum(any(nearest_voxel < 0.5)) || sum(any(nearest_voxel > Dim + 0.5))
    out = 0;
    return
else
    out = mask(convind(nearest_voxel, Dim));
end

end

% if Dim(1) == 1
%     out = out(2:end); %i.e. out(2)
% end
% if any(index > Dim) || any(index < Dim)
%     out = 0;
%     return
% end
