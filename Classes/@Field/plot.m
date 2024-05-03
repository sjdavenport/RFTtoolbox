function out = plot( varargin )
% PLOT( field ) redefines the basic plot routine for scalar fields.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  field   an object of class Field with fiberD = 1 
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
D = varargin{1}.D;


%% Check mandatory input
%--------------------------------------------------------------------------

% Check the field input
if( D ~= 1 )
    error( "The field input must be a field over a line." )
end

% Check the field input
if( varargin{1}.fiberD ~= 1 )
    error( "The field input must be a scalar field, i.e. fiberD = 1." )
end

%% Check optional input
%--------------------------------------------------------------------------

%% Main function
%--------------------------------------------------------------------------
out = plot( varargin{1}.xvals{1}, squeeze( varargin{1}.field ),...
            varargin{2:end} );

return