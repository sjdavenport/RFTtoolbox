function plotPixelSides(mask, linewidth)
% PLOTPIXELSIDES plots lines for each side of each pixel in a 2D binary mask.
%
%   PLOTPIXELSIDES(mask) takes a 2D binary mask as input and displays lines
%   for each side of each filled pixel in the mask. The mask should be a
%   logical or numeric array where 1 represents filled pixels and 0
%   represents empty pixels. The function plots the lines using a blue
%   color. The x-axis represents the column indices, and the y-axis
%   represents the row indices of the mask.
%--------------------------------------------------------------------------
% % Example:
%  mask = [ 0 0 0 0 0;
%           0 0 0 1 0;
%           0 1 0 1 0;
%           0 1 1 1 0;
%           0 0 0 0 0;];
%       plotPixelSides(mask);
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
%     mask = flipud(mask);

if ~exist('linewidth', 'var')
    linewidth = 2;
end
    
    % Get the size of the mask
    [rows, cols] = size(mask);
    
    % Initialize figure and axes
    hold on;
    axis equal;
    
    % Loop over each pixel in the mask
    for i = 1:rows
        for j = 1:cols
            % Check if the pixel is set to 1
            if mask(i, j) == 1
                % Calculate the coordinates of the pixel corners
                x = [j-0.5, j+0.5, j+0.5, j-0.5, j-0.5];
                y = [i-0.5, i-0.5, i+0.5, i+0.5, i-0.5];
                
                % Plot lines for each side of the pixel
                plot(x, y, 'black', 'linewidth', linewidth);
            end
        end
    end
    
    % Set plot limits
    xlim([0.5, cols+0.5]);
    ylim([0.5, rows+0.5]);
%     axis off
    
    % Turn off the hold
    hold off;
end
