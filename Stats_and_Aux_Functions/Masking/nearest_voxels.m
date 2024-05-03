function [ nearby_voxels, nearest_ones ]  = nearest_voxels( point )
% nearestvoxels(point) finds the nearest voxel or voxels to a given point
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  point: a D x 1 vector where D is the number of dimensions, representing
%  the point for which to find the nearest voxels.
%--------------------------------------------------------------------------
% OUTPUT
% nearby_voxels: an D x K matrix, where D is the number of dimensions and
%                K is the number of nearby voxels. Each column contains 
%                the coordinates of a nearby voxel.
% nearest_ones: an D x L matrix, where D is the number of dimensions and 
%               L is the number of nearest voxels. Each column contains the 
%               coordinates of a nearest voxel.
%--------------------------------------------------------------------------
% EXAMPLES
% nearest_voxels([2.4;2])
% nearest_voxels([2.4;1.3])
% [nv, no]= nearest_voxels([2.4;1.3;-0.7])
% [nv, no] = nearest_voxels([2.5;1.2])
% [nv, no] = nearest_voxels([2.5;1.5])
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

nearby_voxels = nearestnonuniquevoxels( point );

nearby_voxels = unique(nearby_voxels','rows')';

D = size(point,1);
if D > 1
    % Find the norm at each voxel using vecnorm
    normdiff = vecnorm(nearby_voxels - point);
else
    % In 1D vecnorm doesn't work and you just need the absolute value
    normdiff = abs(nearby_voxels - point);
end
min_norm = min(normdiff);

nearest_ones = nearby_voxels(:,normdiff < min_norm + 10^(-10));

% % if D == 1
%    nearby_voxels(1) = floor(point);
%    nearby_voxels(2) = ceil(point);
%    nearby_voxels = unique(nearby_voxels);
% elseif D == 2
%    nearby_voxels(1,1:2) = floor(point(1));
%    nearby_voxels(1,3:4) = ceil(point(1));
%    nearby_voxels(2,2:3) = floor(point(2));
%    nearby_voxels(2,[1,4]) = ceil(point(2));
% elseif D == 3
%    nearby_voxels(1,1:4) = floor(point(1));
%    nearby_voxels(1,5:8) = ceil(point(1));
% end


end

function nearby_voxels = nearestnonuniquevoxels( point )
D = size(point,1);
nearby_voxels = zeros(D, 2^(D));

nearby_voxels(1,1:2^(D-1)) = floor(point(1));
nearby_voxels(1,(2^(D-1)+1):2^(D)) = ceil(point(1));

if D > 1
    nearby_voxels(2:end,1:2^(D-1)) = nearestnonuniquevoxels(point(2:end));
    nearby_voxels(2:end,(2^(D-1)+1):2^(D)) = nearby_voxels(2:end,1:2^(D-1));
end

end
%
% tol = 10^(-10);
%
% modpoint5 = mod(point,0.5);
% modpoint5lessthantol = find(modpoint5 < tol);
% mod1 = mod(point, 1);
% mod1lessthantol = find(mod1 < tol);
%
% point5locs = modpoint5lessthantol.*(~mod1lessthantol);
%
% % If no points have a 0.5 in them
% if ~any(point5locs)
%     point = round(point);
%     return
% end
%
% % Obtain the number of dimensions
% D = size(point,1);
%
% % Otherwise return the surrouding voxels
% inds = find(point5locs);
%
% setofnearbyvoxels = cell(1, 2^length(inds));
%
% for d = 1:D
%
% end
%
