function out = applyinterfield( tval, lat_data, xvals_vecs )
% The applyinterfield function is used to interpolate a given field (stored 
% in the lattice field lat_data) at specified values (stored in the tval 
% variable). The function is designed to work in dimensions 1, 2 and 3 and 
% uses the xvals_vecs cell array to determine the positions of the 
% voxels in the lattice field. It uses the nearby_voxels function to select
% the closest voxels to the point of interpolation and applies a linear 
% interpolation kernel to calculate the value of the field at the specified
% point. The output of the function is a single value, corresponding to the 
% interpolated value of the field at the specified point.
%--------------------------------------------------------------------------
% ARGUMENTS
% tval      the t values (an ndim = D by nvalues matrix) at which to evaluate 
%           the field.
% lat_data  the lattice field an array of size Dim.
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
% 
%--------------------------------------------------------------------------
% EXAMPLES
% lat_data = [0,100;1,3];
% tval = [1,1]';
% out = applyinterfield( tval, lat_data, {1:2,1:2} )
%--------------------------------------------------------------------------
% Samuel Davenport
%--------------------------------------------------------------------------
D = length(xvals_vecs);
increm = zeros(D,1);
for d = 1:D
   increm(d) = xvals_vecs{d}(2) - xvals_vecs{d}(1);
end

%%  Main Function Loop
%--------------------------------------------------------------------------
pointdividedbyincrem = tval./increm;
nearbyvoxels = nearest_voxels( pointdividedbyincrem );

if size(nearbyvoxels, 2) == 1
    out = lat_data(xvals_vecs{1}(nearbyvoxels(1)), ...
        xvals_vecs{2}(nearbyvoxels(2)));
    return
end

if D == 1
    x = xvals_vecs{1}(max(nearbyvoxels(1,:))) - tval(1);
    x = 1-x;
    V1 = lat_data(nearbyvoxels(1,1));
    V2 = lat_data(nearbyvoxels(1,2));
    out = (1-x)*V1 + x*V2;
end

if D == 2
    x = xvals_vecs{1}(max(nearbyvoxels(1,:))) - tval(1);
    y = xvals_vecs{2}(max(nearbyvoxels(2,:))) - tval(2);
    x = 1-x;
    y = 1-y;
    
    V1 = lat_data(nearbyvoxels(1,1), nearbyvoxels(2,1));
    V2 = lat_data(nearbyvoxels(1,2), nearbyvoxels(2,2));
    V3 = lat_data(nearbyvoxels(1,3), nearbyvoxels(2,3));
    V4 = lat_data(nearbyvoxels(1,4), nearbyvoxels(2,4));
    out = (1-x)*(1-y)*V1 + x*(1-y)*V2 + (1-x)*y*V3 + x*y*V4;
end

if D == 3
    x = xvals_vecs{1}(max(nearbyvoxels(1,:))) - tval(1);
    y = xvals_vecs{2}(max(nearbyvoxels(2,:))) - tval(2);
    z = xvals_vecs{3}(max(nearbyvoxels(3,:))) - tval(3);
    x = 1-x;
    y = 1-y;
    z = 1-z;
    
    V1 = lat_data(nearbyvoxels(1,1), nearbyvoxels(2,1), nearbyvoxels(3,1));
    V2 = lat_data(nearbyvoxels(1,2), nearbyvoxels(2,2), nearbyvoxels(3,2));
    V3 = lat_data(nearbyvoxels(1,3), nearbyvoxels(2,3), nearbyvoxels(3,3));
    V4 = lat_data(nearbyvoxels(1,4), nearbyvoxels(2,4), nearbyvoxels(3,4));
    V5 = lat_data(nearbyvoxels(1,5), nearbyvoxels(2,5), nearbyvoxels(3,5));
    V6 = lat_data(nearbyvoxels(1,6), nearbyvoxels(2,6), nearbyvoxels(3,6));
    V7 = lat_data(nearbyvoxels(1,7), nearbyvoxels(2,7), nearbyvoxels(3,7));
    V8 = lat_data(nearbyvoxels(1,8), nearbyvoxels(2,8), nearbyvoxels(3,8));
    out = (1-x)*(1-y)*(1-z)*V1 + x*(1-y)*(1-z)*V2 + (1-x)*y*(1-z)*V3 + x*y*(1-z)*V4 + ...
          (1-x)*(1-y)*z*V5 + x*(1-y)*z*V6 + (1-x)*y*z*V7 + x*y*z*V8;
end

end

