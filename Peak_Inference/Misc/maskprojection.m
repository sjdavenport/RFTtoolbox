function projected_point = maskprojection( point, mask )
% MASKPROJECTION( point, mask ) takes a point near (within half a voxel of) 
% the boundary of a mask and projects it to the mask.
%--------------------------------------------------------------------------
% ARGUMENTS
%
%--------------------------------------------------------------------------
% OUTPUT
%
%--------------------------------------------------------------------------
% EXAMPLES
% %1D
% maskprojection( 0.1, ones(1,4))
% maskprojection( 4.75, ones(1,4))
%
% %2D
% mask = [0,1,1;0,0,1];
% maskprojection( [1.4,1.4]', mask )
%
% mask = [0,0,1;0,0,1;0,0,1];
% maskprojection([2.49,2.49]', mask)
%
% maskprojection([10.7,1.7]',ones(10,10))
% maskprojection([9.7,10.51]', ones(10,10))
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

if inmask(point, mask)
    projected_point = point;
    return
end

Dim = size(mask);
D = length(Dim);
if D == 2 && length(point) == 1
    D = 1;
end
if size(point, 1) ~= D
    error('the point must be a column vector with the same dimension as the mask')
end
point_floor = floor(point);
point_ceil = ceil(point);

surrounding_valsvecs = cell(1,D);
for d = 1:D
    surrounding_valsvecs{d} = [point_floor(d),point_ceil(d)];
end

surrounding_voxels = xvals2voxels(surrounding_valsvecs);
mask_voxels_ind = zeros(1, size(surrounding_voxels,2));
for I = 1:length(surrounding_voxels)
    mask_voxels_ind(I) = inmask(surrounding_voxels(:,I), mask);
end
if sum(mask_voxels_ind) == 0
    error('the point must lie within 1/2 a voxel of the mask!')
end
mask_voxels_ind = logical(mask_voxels_ind); %So that you can subset the matrix!

surrounding_voxels_inmask = surrounding_voxels(:, mask_voxels_ind);
[~, nearest_voxel_inmask_ind] = min(sum((surrounding_voxels_inmask - point).^2,1));
nearest_voxel_inmask = surrounding_voxels_inmask(:, nearest_voxel_inmask_ind);
points_on_boundary_of_nearest_voxel_valsvecs = cell(1,D);
for d = 1:D
    points_on_boundary_of_nearest_voxel_valsvecs{d} = [nearest_voxel_inmask(d) - 0.5, nearest_voxel_inmask(d) + 0.5];
end
points_on_boundary_of_nearest_voxel = xvals2voxels(points_on_boundary_of_nearest_voxel_valsvecs);
if D == 1
    [~, sorted_ind] = sort((points_on_boundary_of_nearest_voxel - point).^2);
    projected_point = points_on_boundary_of_nearest_voxel(sorted_ind(1));
elseif D == 2
    [~, sorted_ind] = sort(sum((points_on_boundary_of_nearest_voxel - point).^2));
    nearest_corner = points_on_boundary_of_nearest_voxel(:, sorted_ind(1));
    second_nearest_corner = points_on_boundary_of_nearest_voxel(:, sorted_ind(2));
    if nearest_corner(1) == second_nearest_corner(1)
        projected_point = [ nearest_corner(1), point(2)]';
        if nearest_corner(1) < point(1) && ~inmask(projected_point-[0.00001,0]', mask)
            projected_point = nearest_corner;
        elseif nearest_corner(1) > point(1) && ~inmask(projected_point+[0.00001,0]', mask)
            projected_point = nearest_corner;
        end
    else
        projected_point = [ point(1), nearest_corner(2)]';
        if nearest_corner(2) < point(2) && ~inmask(projected_point-[0,0.00001]', mask)
            projected_point = nearest_corner;
        elseif nearest_corner(2) > point(2) && ~inmask(projected_point+[0,0.00001]', mask)
            projected_point = nearest_corner;
        end
    end
elseif D == 3
%     [~, sorted_ind] = sort(sum((points_on_boundary_of_nearest_voxel - point).^2));
    %             nearest_corners = points_on_boundary_of_nearest_voxel(:, sorted_ind(1:4));
    %NEed to work This Bit oout1!
end

end

% function out = inmask(point, mask)
%     Dim = size(mask);
%     nearest_voxel = round(point);
%     if sum(any(nearest_voxel < 0)) || sum(any(nearest_voxel > Dim))
%         out = 0;
%     else
%         out = mask(convind(nearest_voxel, Dim));
%     end
% end
