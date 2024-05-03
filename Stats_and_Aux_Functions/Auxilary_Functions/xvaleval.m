function evaluated_points = xvaleval( points, xvals )
% XVALEVALS( points, xvals ) takes the point allocations on a lattice and
% finds their values in terms of the xvals coordinates
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  points   a D by npoints matrix
%  xvals    a D-dimensional cell array whose entries are vectors giving 
%           the xvalues at each each dimension along the lattice. It 
%               assumes a regular, rectangular lattice (though within a given
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
%  evaluated_points    a D by npoints matrix where each column is a point
%                      that has been evaluated using the givne xvals
%--------------------------------------------------------------------------
% EXAMPLES
% %% 1D
% xvaleval(6, 1:0.5:10)
% %% 2D
% xvaleval([2;3], 1:0.5:3)
% xvaleval([2;3], 1:0.5:3)
% xvaleval([2;3], {1:0.5:3,1:0.5:3})
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Obtain the number of dimensions
D = size(points,1);

% Convert xvals to a cell if it is not.
if isnumeric(xvals)
    xvals = {xvals};
end

% Ensure that xvals has the right dimension
if length(xvals) == 1
    xvals = repmat(xvals, 1, D);
elseif length(xvals) < D
    error('xvals and points must have the same number of dimensions')
end

%%  Main Function Loop
%--------------------------------------------------------------------------
evaluated_points = zeros(size(points));

for d = 1:D
    dthpointdim = points(d, :);
    if max(dthpointdim) > length(xvals{d})
        error('The points must be inside xvals')
    end
    evaluated_points(d, :) = xvals{d}(dthpointdim);
end

end

