function out = imagesc( field, slice, subj )
% PLOT( field ) redefines the basic plot routine for scalar fields.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  field   an object of class Field with fiberD = 1
% Optional
%  slice   an 1 x field.D vector containing the location of the voxel at
%          which the 2D slices should be plotted. Note that two components
%          need to be NaN. These components define the plotted slice.
%  subj    an integer indicating which subject should be plotted. Default
%          1.
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

% Check the field input
if( field.fiberD ~= 1 )
    error( "The field input must be a scalar field, i.e. fiberD = 1." )
end

if( D < 2 )
    error( "The field input must have dimension >= 2." )
end

%% Check optional input
%--------------------------------------------------------------------------

% Check slice input
if ~exist( 'slice', 'var' )
    slice    = ceil( smask / 2 );
    slice(1) = NaN;
    slice(2) = NaN;
elseif isnumeric( slice )
    % Check the field input
    if sum( isnan( slice ) ) ~= 2 && length( slice ) ~= D
        error( "You need to include exactly two NaNs in the 1 x field.D slice vector." )
    end
else
    error( "slice needs to be a numeric vector" )
end

% Check slice input
if ~exist( 'subj', 'var' )
    subj = 1;
end

%% Main function
%--------------------------------------------------------------------------

% create index for the chosen slice
index = cell( [ 1 D ] );
ind_NaN = NaN * [ 1 1 ];
count = 1;
for d = 1:D
    if ~isnan( slice(d) )
        index{d} = slice(d);
    else
        index{d} = ':';
        ind_NaN(count)  = d;
        count = 2;
    end
end

xvals1 = field.xvals{ ind_NaN(1) };
xvals2 = field.xvals{ ind_NaN(2) };

% number of voxels between shown xvals
n = 9;

% plot the image
out = imagesc( squeeze( field.field( index{:}, subj ) ) );
set(gca, 'XTick', 1:n:smask( ind_NaN(1) ),...
         'XTickLabel', xvals1(1:n:smask( ind_NaN(1) ) ) )
set(gca, 'YTick', 1:n:smask( ind_NaN(2) ),...
         'YTickLabel', xvals2(1:n:smask( ind_NaN(2) ) ) )

return