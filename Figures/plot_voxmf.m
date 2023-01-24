function [ ] = plot_voxmf( field )
% NEWFUN( var1, var2, opt1 ) is described here [...]
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   var1   this is a mandatory variable.
%   var2   this is a mandatory variable.
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
%%%%%% HEADLINE 
%%% SUBHEADLINE
% Line or short paragraph comment

%% Add/check optional values
%--------------------------------------------------------------------------

% This kind of code with exists is better than using nargin < xy, since
% Parameter can be easily permuted
if ~exist( 'opt1', 'var' )
   % Default option of opt1
   opt1 = 0;
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
    plotcube(delta, mpoints(i,:,:),.8,[1 0 0]); hold on
end
xlim([min(X(:)), max(X(:))])
ylim([min(Y(:)), max(Y(:))])
zlim([min(Z(:)), max(Z(:))])
axis square
return