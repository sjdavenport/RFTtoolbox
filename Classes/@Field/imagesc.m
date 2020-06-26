function out = imagesc( field, slice, ntics )
% imagesc( field, slice, subj, ntics ) redefines the imagesc function for
% visualization of objects of class Field.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  field   an object of class Field with fiberD = 1
% Optional
%  slice   an 1 x (field.D + fiberD) cell array containing the location of
%          the voxel at which the 2D slices should be plotted.
%          Note that two components need to be NaN. These components define
%          the spatial part of the plotted slice.
%  ntics   an integer denoting the number of tics shown in x and y direction.
%--------------------------------------------------------------------------
% OUTPUT
%   out  an handle of the imagesc plot
% -------------------------------------------------------------------------
% DEVELOPER TODOs:
%--------------------------------------------------------------------------
% EXAMPLES
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow
%--------------------------------------------------------------------------

%% Get constants
%--------------------------------------------------------------------------

% Dimension of the domain
D = field.D;

% Get number of subjects/samples
smask = field.masksize;

%% Check mandatory input
%--------------------------------------------------------------------------

if( D < 2 )
    error( "The field input must have dimension >= 2." )
end

%% Check optional input
%--------------------------------------------------------------------------

% Check slice input
if ~exist( 'slice', 'var' )
    slice    = num2cell( [ ceil( smask / 2 ), ones( [1, field.fiberD ] ) ] );
    slice{1} = NaN;
    slice{2} = NaN;
elseif iscell( slice )
    % Turn into vector
    test = [ slice{:} ];
    test = test( 1:D );
    
    % Check the field input
    if sum( isnan( test ) ) ~= 2 && length( test ) ~= D
        error( "You need to include exactly two NaNs in the 1 x field.D slice vector." )
    end
else
    error( "slice needs to be a numeric vector" )
end

% Check slice input
if ~exist( 'ntics', 'var' )
    ntics = 10;
end

%% Main function
%--------------------------------------------------------------------------

% Create index for the chosen slice
count = 1;
ind = [ NaN NaN ];
for d = 1:D
    if count == 1 && isnan( slice{d} )
        end1 = field.masksize(d);
        slice{d} = end1:-1:1;
        % Reverse x dimension
        xvals1 = field.xvals{ d };
        % Get voxelsize
        dx1    = xvals1(2) - xvals1(1);
        count  = 2;
        ind(1) = d;
    elseif count == 2 && isnan( slice{d} )
        end2 = field.masksize(d);
        slice{d} = ':';
        xvals2 = field.xvals{ d };
        % Get voxelsize
        dx2    = xvals2(2) - xvals2(1);
        count  = 3;
        ind(2) = d;
    end
end

% Reduce input field to scalar field with 2D domain
field = Subfield( field, slice );
slice = slice(ind);


% number of voxels between shown xvals
nx    = ceil( end1 / ( ntics) );
xtics = ( nx ) : nx : end1;
ny    = ceil( end2 / ( ntics) );
ytics = ( ny ) : ny : end2;

% plot the image
out = imagesc( squeeze( field.field( slice{:} ) ) );
set(gca, 'YTick', xtics ,...
         'YTickLabel', xvals1( xtics( end:-1:1 ) ) )
set(gca, 'XTick', ytics ,...
         'XTickLabel', xvals2( ytics ) )

daspect( [ dx1 dx2 max(dx1, dx2) ])     

return