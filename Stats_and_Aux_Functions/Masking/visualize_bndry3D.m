function [ ] = visualize_bndry3D( mask, resadd, types, pts_size, angle )
% visualize_bndry( mask, resadd, types, pts_size, angle ) visualizes the
% resolution increased boundary from mask_highres.m as well as its subparts,
% which are computed using bndry_voxel.m.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   mask    a logical T_1 x ... x T_D array.
%
% Optional
%   resadd  the amount of equidistant voxels introduced inbetween the
%           voxels  
%   types   a string or vector of strings indicating which part of the
%           boundary should be obtained.
%           For D = 2 further options are
%               - 'y', which returns all boundary segments with fixed
%                      x-value
%               - 'x', which returns all boundary segments with fixed
%                      y-value
%           For D = 3 further options are
%               - 'xy', which returns all voxels of boundary faces with
%                       fixed z-value
%               - 'xz', which returns all voxels of boundary faces with
%                       fixed y-value
%               - 'yz', which returns all voxels of boundary faces with
%                       fixed x-value
%               - 'x',  which returns all edges in x-direction
%               - 'y',  which returns all edges in y-direction
%               - 'z',  which returns all edges in z-direction
%   pts_size a numeric to increase or decrease the size of the points.
%            Default 50.
%   angle    a 1 x 2 vector specifying the 3D view. First component rotates
%            around the vertical z-axis, second component around an axis in
%            the x-y plane. Default [ 42 20 ].
%
%--------------------------------------------------------------------------
% OUTPUT
%   Outputs an image showing the resolution increased boundary
%--------------------------------------------------------------------------
% DEVELOPER TODOs: - add 2D option
%--------------------------------------------------------------------------
% % EXAMPLES
% %% % Show box example all boundaries
% visualize_bndry3D( true([4 4 4]), 1, ["yz", "xz", "xy"], 50 )
% title("All boundary voxels colored by constant coordinate")
% %% % Show box example only one boundary
% visualize_bndry3D( true([4 4 4]), 1, ["yz", "x"], 50 )
% title("y-z-plane faces of the boundary and x-edges")
% %% %% Spherical object
% % Create a mask and show it
% siz = 3;
% dx  = 0.5;
% [x,y,z] = meshgrid( -siz:dx:siz, -siz:dx:siz, -siz:dx:siz );
% xvals = [x(:), y(:), z(:)]';
% h     = reshape( GkerMV( xvals, 5 ), size(x) );
% mask  = logical( h > 0.003 );
% clear h
% %% % resadd = 1
% visualize_bndry3D( mask, 1, [ "xy", "yz", "xz" ], 40 )
% title("All boundary voxels colored by constant coordinate")
% 
% %% % resadd = 3
% visualize_bndry3D( mask, 3, [ "x", "y", "z" ], 40 )
% title("All edges colored by coordinate")
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow  
%--------------------------------------------------------------------------


%% Check mandatory input and get important constants
%--------------------------------------------------------------------------

if ~length( size( mask ) ) == 3
    error( "mask must be 3D array." );
end

%% Add/check optional values
%--------------------------------------------------------------------------

% This kind of code with exists is better than using nargin < xy, since
if ~exist( 'resadd', 'var' )
   % Default option of opt1
   resadd = 1;
end

% Get distance between high res points
dx = 1 / ( resadd + 1 );

% Default plot all boundary values
if ~exist( 'types', 'var' )
   % Default option of opt1
   types = ["xy", "xz", "yz"];
end

% Default plot all boundary values
if ~exist( 'angle', 'var' )
   % Default option of opt1
   angle = [ 42 20 ];
end

% Default plot all boundary values
if ~exist( 'pts_size', 'var' )
   % Default option of opt1
   pts_size = 50;
end

%% Main function  
%--------------------------------------------------------------------------

% Define a color structure
colors = [ "b", "y", "k" ];
spts = pts_size * [ 3, 1, 0.5 ];
           
% Get the highres mask
mask_hr = mask_highres( mask, resadd );

% Get the coordinates of the original voxels
[m1, m2, m3] = ind2sub( size( mask ), find( mask == 1 ) );
mask_pts = [m1, m2, m3];
mask_voxsize = [1, 1, 1];

% Get the chosen boundary type
[ bndry, ~ ] = bndry_voxels( logical( mask_hr ), types );

if length( types ) == 1
    % bdry voxel coordinates
    [m1, m2, m3] = ind2sub( size( mask_hr ), find( bndry == 1 ) );

    figure(1), clf, hold on
    % Plot the original voxels
    voxel_image( mask_pts, mask_voxsize, "red", 0.65 );
    % Plot the bdry locations
    scatter3( m1*dx +(0.5-dx), m2*dx+(0.5-dx), m3*dx+(0.5-dx),...
              pts_size, 'blue', 'filled');
    
else
    % Initialize the figure
    figure(1), clf, hold on
    % Plot the boundary points for the chosen faces
    k = 0;
    for type = types
        k = k+1;
        % bdry voxel coordinates
        [m1, m2, m3] = ind2sub( size( mask_hr ), find( bndry.(type) == 1 ) );
        % Plot the original voxels
        voxel_image( mask_pts, mask_voxsize, "red", 0.65 );
        % Plot the bdry locations
%         scatter3( m1*dx +(0.5-dx), m2*dx+(0.5-dx), m3*dx+(0.5-dx),...
%                   spts(k), colors(k), 'filled' );
        scatter3( m1*dx +(0.5-dx), m2*dx+(0.5-dx), m3*dx+(0.5-dx),...
                  spts(k), 'yellow', 'filled' );
    end
end
% Change the view angle
view( angle )
axis equal

return