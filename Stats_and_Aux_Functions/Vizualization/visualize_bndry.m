function visualize_bndry( mask, resadd, types, pts_size, center_vox_size_ratio, alpha_level, angle, linewidth )
% visualize_bndry( mask, resadd, types, pts_size, angle ) visualizes the
% resolution increased boundary from mask_highres.m as well as its subparts,
% which are computed using bndry_voxel.m.
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
%--------------------------------------------------------------------------
% OUTPUT
%   Outputs an image showing the resolution increased boundary
%--------------------------------------------------------------------------
% DEVELOPER TODOs: - add 2D option
%--------------------------------------------------------------------------
% % EXAMPLES
% %% % Show box example all boundaries
% visualize_bndry( true([4 4 4]), 1, ["yz", "xz", "xy"], 50 )
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
% visualize_bndry( mask, 1, [ "xy", "yz", "xz" ], 40 )
% title("All boundary voxels colored by constant coordinate")
%
% %% % resadd = 3
% visualize_bndry( mask, 3, [ "x", "y", "z" ], 40 )
% title("All edges colored by coordinate")
%
% % % Double Convex Edge
% mask = false([3 3 3]);
% mask( 2, 2, 2 ) = 1;
% mask( 2, 3, 3 ) = 1;
% types = ["x"];
% pts_size = 100;
% resadd = 1;
% visualize_bndry( mask, resadd, types, pts_size)
% view( [ 30 20 ] )
%--------------------------------------------------------------------------
% AUTHORS: Fabian Telschow and Samuel Davenport
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

if ~exist('linewidth', 'var')
    linewidth = 1;
end

if ~exist( 'alpha_level', 'var' )
    % Default option of opt1
    alpha_level = 0.65;
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

if ~exist('center_vox_size_ratio', 'var')
    center_vox_size_ratio = 4;
end

%% Main function
%--------------------------------------------------------------------------

% Define a color structure
colors = [ "b", "y", "k" ];
% spts = pts_size * [ 3, 1, 0.5 ];
spts = repmat(pts_size, 1, 3);

% Get the highres mask
% A temporary fix for the [1,1,1] case which doesn't work.
D = length(size(mask));
if isequal(size(mask), [1,1])
    mask_hr = ones(repmat(resadd + 2, 1,3));
    D = 3;
else
    mask_hr = mask_highres( mask, resadd );
end

% Get the coordinates of the original voxels
if D == 3
    [a, b, c] = ind2sub( size( mask ), find( mask == 1 ) );
    mask_pts = [a, b, c];
    mask_voxsize = [1, 1, 1];
elseif D == 2
    [a, b] = ind2sub( size( mask ), find( mask == 1 ) );
    mask_pts = [a, b];
    mask_voxsize = [1, 1];
    center_vox_size_ratio = 1;
end

% Get the chosen boundary type
if ~strcmp(types, "all")
    [ bndry, ~ ] = bndry_voxels( logical( mask_hr ), types );
end

if length( types ) == 1
    % bdry voxel coordinates
    
    if strcmp(types, "all")
        a = spacep(a, resadd) + (resadd + 1)/2;
        b = spacep(b, resadd) + (resadd + 1)/2;
        if D == 3
            [m1, m2, m3] = ind2sub( size( mask_hr ), find( mask_hr == 1 ) );
            c = spacep(c, resadd) + (resadd + 1)/2;
        else
            [m1, m2] = ind2sub( size( mask_hr ), find( mask_hr == 1 ) );
        end
    else
        if D == 3
            [m1, m2, m3] = ind2sub( size( mask_hr ), find( bndry == 1 ) );
        else
            [m1, m2] = ind2sub( size( mask_hr ), find( bndry == 1 ) );
        end
    end
    
    figure(1), clf, hold on
    % Plot the original voxels
    if D == 3
        voxel_image( mask_pts, mask_voxsize, "red", alpha_level );
    else
        voxel_image_2D( flipud(rot90(mask,1)),  alpha_level, linewidth, "red");
    end
    % Plot the bdry locations
    hold on
    if D == 3
        if resadd > 0
            scatter3( m1*dx +(0.5-dx), m2*dx+(0.5-dx), m3*dx+(0.5-dx),...
                pts_size, 'yellow', 'filled');
        end
        
        if strcmp(types, "all")
            scatter3( a*dx +(0.5-dx), b*dx+(0.5-dx), c*dx+(0.5-dx),...
                pts_size*center_vox_size_ratio, 'yellow', 'filled');
        end
    elseif D == 2
        if resadd > 0
            scatter( m1*dx +(0.5-dx), m2*dx+(0.5-dx), pts_size, 'yellow', 'filled');
        end
        if strcmp(types, "all")
            scatter( a*dx +(0.5-dx), b*dx+(0.5-dx), pts_size*center_vox_size_ratio, 'yellow', 'filled');
        end
    end
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
        voxel_image( mask_pts, mask_voxsize, "red", alpha_level );
        % Plot the bdry locations
        %         scatter3( m1*dx +(0.5-dx), m2*dx+(0.5-dx), m3*dx+(0.5-dx),...
        %                   spts(k), colors(k), 'filled' );
        scatter3( m1*dx +(0.5-dx), m2*dx+(0.5-dx), m3*dx+(0.5-dx),...
            spts(k), 'yellow', 'filled' );
        %         pause
    end
end
% Change the view angle
if D == 3
view( angle )
end
axis equal
axis off

return