function out = plot( field, slice )
% PLOT( field ) redefines the basic plot routine for scalar fields.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  field   an object of class Field with fiberD = 1
% Optional
%  slice   an 1 x field.D vector containing the location of the voxel at
%          which the 1D slices should be plotted. Note that one component
%          needs to be NaN. This component is considered to be sliced.     
%--------------------------------------------------------------------------
% OUTPUT
%   out  an handle of the plot
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

%% Check optional input
%--------------------------------------------------------------------------

% Check slice input
if ~exist( 'slice', 'var' )
    slice    = ceil( smask / 2 );
    slice(1) = NaN;
elseif isnumeric( slice )
    % Check the field input
    if sum( isnan( slice ) ) ~= 1 && length( slice ) ~= D
        error( "You need to include exactly one NaN in the 1 x field.D slice vector." )
    end
else
    error( "slice needs to be a numeric vector" )
end


%% Main function
%--------------------------------------------------------------------------

% create index for the chosen slice
index = cell( [ 1 D ] );
for d = 1:D
    if ~isnan( slice(d) )
        index{d} = slice(d);
    else
        index{d} = ':';
        xvals  = field.xvals{d};
    end
end

out = plot( xvals, squeeze( field.field( index{:}, : ) ) );

return