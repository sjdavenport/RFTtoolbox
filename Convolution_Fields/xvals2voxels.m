function xvaluesatvoxels = xvals2voxels( xvals_vecs, D )
% XVALS2VOXELS transforms the vectors giving the coordinates of the sides
% of the lattice into values of all of the points within the lattice.
%--------------------------------------------------------------------------
% ARGUMENTS
% xvals_vecs    a D-dimensional cell array whose entries are vectors giving the
%               xvalues at each each dimension along the lattice. It assumes
%               a regular, rectangular lattice (though within a given
%               dimension the voxels can be spaced irregularly).
%               I.e suppose that your initial lattice grid is a
%               4by5 2D grid with 4 voxels in the x direction and 5 in
%               the y direction. And that the x-values take the values:
%               [1,2,3,4] and the y-values take the values: [0,2,4,6,8].
%               Then you would take xvals_vecs = {[1,2,3,4], [0,2,4,6,8]}.
%               The default is to assume that the spacing between the
%               voxels is 1. If only one xval_vec direction is set the
%               others are taken to range up from 1 with increment given by
%               the set direction.
%--------------------------------------------------------------------------
% OUTPUT
% xvaluesatvoxels
%--------------------------------------------------------------------------
% EXAMPLES
% xvals_vecs = {[1,2,3,4], [0,2,4,6]};
% xvaluesatvoxels = xvals2voxels({[1,2,3,4], [0,2,4,6]})
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

if ~iscell(xvals_vecs)
    xvals_vecs = {xvals_vecs};
end

if nargin < 2
    D = length(xvals_vecs);
else
    xvals_vecs = repmat(xvals_vecs, 1, D);
end

if D == 1
    xvaluesatvoxels = xvals_vecs{1};
elseif D == 2
    [x, y] = ndgrid(xvals_vecs{1},xvals_vecs{2});
    %     [x, y] = ndgrid(1:Ydim(1),1:Ydim(2));
    xvaluesatvoxels = [x(:),y(:)]';
elseif D == 3
    [x, y, z] = ndgrid(xvals_vecs{1},xvals_vecs{2},xvals_vecs{3});
    xvaluesatvoxels = [x(:),y(:),z(:)]';
end

end

