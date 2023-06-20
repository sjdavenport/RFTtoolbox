function [ out ] = voxel_image_2D( mask, alpha_level, linewidth, color )
% NEWFUN
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% % EXAMPLES
% mask = [ 0 0 0 0 0;
%           0 0 0 1 0;
%           0 1 0 1 0;
%           0 1 1 1 0;
%           0 0 0 0 0;];
% voxel_image_2D(mask)
% mask = ones([2,2]);
% voxel_image_2D(mask, 0.5)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'color', 'var' )
   % Default value
   color = 'red';
end

if ~exist( 'alpha', 'var' )
   % Default value
   alpha_level = 0.5;
end


if ~exist( 'linewidth', 'var' )
   % Default value
   linewidth = 1;
end


%%  Main Function Loop
%--------------------------------------------------------------------------
viewdata( mask, ones(size(mask)), {mask, 1-mask}, {color, 'white'}, 1, [], alpha_level)
plotPixelSides(mask, linewidth);

end

