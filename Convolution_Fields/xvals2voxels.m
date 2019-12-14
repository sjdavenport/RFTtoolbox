function xvalues_at_voxels = xvals2voxels( xvals_vecs )
% XVALS2VOXELS transforms the vectors giving the coordinates of the sides
% of the lattice into values of all of the points within the lattice.
%--------------------------------------------------------------------------
% AUTHOR: Samuel J. Davenport

D = length(xvals_vecs);

if D == 1
    xvalues_at_voxels = xvals_vecs{1}';
elseif D == 2
    [x, y] = ndgrid(xvals_vecs{1},xvals_vecs{2});
    %     [x, y] = ndgrid(1:Ydim(1),1:Ydim(2));
    xvalues_at_voxels = [x(:),y(:)];
elseif D == 3
    [x, y, z] = ndgrid(xvals_vecs{1},xvals_vecs{2},xvals_vecs{3});
    xvalues_at_voxels = [x(:),y(:),z(:)];
end

end

