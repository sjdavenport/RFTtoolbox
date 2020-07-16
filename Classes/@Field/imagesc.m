function out = imagesc( field, ntics )
% imagesc( field, slice, subj, ntics ) redefines the imagesc function for
% visualization of objects of class Field.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  field   an object of class Field with fiberD = 1
% Optional
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

if( D ~= 2 )
    error( "The field input must have dimension = 2." )
end

if( field.fiberD ~= 1 )
    error( "The field input must have fiberD = 1." )
end

%% Check optional input
%--------------------------------------------------------------------------

% Check slice input
if ~exist( 'ntics', 'var' )
    ntics = 10;
end

%% Main function
%--------------------------------------------------------------------------

% Create index for the chosen slice
end1 = smask(1);
end2 = smask(2);

% xvals
xvals1 = field.xvals{ 1 };
xvals2 = field.xvals{ 2 };

% Get voxelsize
dx1 = xvals1(2) - xvals1(1);
dx2 = xvals2(2) - xvals2(1);

I = field.field( end1:-1:1, : );

% number of voxels between shown xvals
nx    = ceil( end1 / ( ntics) );
xtics = ( nx ) : nx : end1;
ny    = ceil( end2 / ( ntics) );
ytics = ( ny ) : ny : end2;

% plot the image
out = imagesc( squeeze( I ) );
set(gca, 'YTick', xtics ,...
         'YTickLabel', xvals1( xtics( end:-1:1 ) ) )
set(gca, 'XTick', ytics ,...
         'XTickLabel', xvals2( ytics ) )

daspect( [ dx1 dx2 max(dx1, dx2) ])     

return