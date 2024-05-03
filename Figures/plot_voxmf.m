function [ ] = plot_voxmf( field, col )
% NEWFUN( var1, var2, opt1 ) is described here [...]
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  field   an object of class Field with field.D = 3.
%   col    a 1x3 vector with values between 0 and 1, giving the RGB color
%          values of the box.
%
% Optional
%   opt1   this is an optional parameter. Default 0.
%
%--------------------------------------------------------------------------
% OUTPUT
%   out1  this is the first output
%   out2  this is the second output
%--------------------------------------------------------------------------
% DEVELOPER TODOs:
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow  
%--------------------------------------------------------------------------


%% Check mandatory input and get important constants
%--------------------------------------------------------------------------
if ~isa(field, "Field")
    error( "The input field must be of Class Field." )
end

if field.D ~= 3
    error( "The input field must have dimension 3." )
end
%% Add/check optional values
%--------------------------------------------------------------------------

% This kind of code with exists is better than using nargin < xy, since
% Parameter can be easily permuted
if ~exist( 'col', 'var' )
   % default number of bootstrap replicates
   col = [70/255 0 150/255];
end


%% Main function  
%--------------------------------------------------------------------------

% Get an xvals array
[X,Y,Z] = meshgrid(field.xvals{1}, field.xvals{2}, field.xvals{3});
Grid = [X(:), Y(:), Z(:)];

% Get mask as a vector
mask_vec = field.mask(:);

% Middle points of the boxes of the voxel manifold
mpoints = Grid(mask_vec,:,:);

% Get the lengths of the voxel manifold boxes
delta = NaN * ones([1, 3]);
for d = 1:3
    delta(d) =min(diff(field.xvals{d}));
end

for i = 1:size(mpoints,1)
    loader(i, size(mpoints,1), 'Percent of the rendering completed:');
    plotcube(delta, mpoints(i,:,:), .7, col); hold on
end
xlim([min(X(:)), max(X(:))])
ylim([min(Y(:)), max(Y(:))])
zlim([min(Z(:)), max(Z(:))])
axis square
grid on
return