function xvals_vecs = getxvalsvecs( Dim, resadd, enlarge, adjust_kernel )
% GETXVALSVECS( Dim, resadd, enlarge ) obtains the xvals_vecs from the
% dimensions of the data and the amount of resolution
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   Dim       a length D vector giving the dimensions of the data
%   resadd    the resolution
% 
% Optional 
%   enlarge   the number of voxels with which to enlarge the mask. 
%             Default is 0
%   adjust_kernel   the amount by which to shift the kernel if necessary.
%                   Default is 0
%--------------------------------------------------------------------------
% OUTPUT
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
% EXAMPLES
% %% 1D
% xvals_vecs = getxvalsvecs( 5, 1, 1 ); xvals_vecs{1};
%
% %% 2D
% xvals_vecs = getxvalsvecs( [5,10], 1 )
% xvals_vecs{1}, xvals_vecs{2}
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

% Default the enlargement to be 0
if ~exist('enlarge', 'var')
    enlarge = 0;
end

% get the number of dimensions
D = length(Dim);

% Default the adjust_kernel to be 0
if ~exist('adjust_kernel', 'var')
    adjust_kernel = zeros(D,1);
end

% get the difference between voxels with resolution increase
dx = 1 / ( resadd + 1 );

%% main loop
%--------------------------------------------------------------------------
xvals_vecs  = cell( 1, D );
for d = 1:D
    xvals_vecs{d} = ( ( 1 - enlarge*dx ):dx:( Dim(d) + enlarge*dx ) )...
                                                        + adjust_kernel(d);
end

end

