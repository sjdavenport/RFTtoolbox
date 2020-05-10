function point_out = spacep( point, spacing )
% SPACEP( POINT, SPACING) records the values of a the location of points in
% a convolution field that has been evaluated on a lattice with given
% spacing
%--------------------------------------------------------------------------
% ARGUMENTS
% point     a D (number of dimensions) by nvals matrix of points
% spacing   the spacing usd between the voxels
%--------------------------------------------------------------------------
% OUTPUT
% point_out  the indices of the converted points
%--------------------------------------------------------------------------
% EXAMPLES
% spacep( [3,3,3]', 0.05)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
inverse_spacing = floor(1/spacing);
mod_spacing = 1/inverse_spacing;

if ~(norm(round( point/mod_spacing ) -  point/mod_spacing) < 10^(-5))
    error('This point is not contained in the spaced lattice')
end

point_out = (point - 1)*inverse_spacing + 1;

end

