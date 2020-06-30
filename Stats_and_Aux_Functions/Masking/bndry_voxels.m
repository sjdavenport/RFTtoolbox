function [ bndry, weights, angle, orientation ] = bndry_voxels( mask, type )
% This function computes boundary voxels and weights for a given mask. Note
% that the weights only make sense for masks which are resolution increased
% by an odd number using mask_highres.m.
% The weights are used in LKC estimation and assume trapozoidal
% integration.
% The boundary is splitted into several subparts representing the faces and
% edges for specific fixed coordinates.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   mask     a logical T_1 x ... x T_D array.
% Optional
%   type  a string or vector of strings indicating which part of the
%            boundary should be obtained.
%            For arbitrary D 'full', which is the default value returns all
%            boundary voxels.
%            For D = 2 further options are
%               - 'y', which returns all boundary segments with fixed
%                      x-value
%               - 'x', which returns all boundary segments with fixed
%                      y-value
%            For D = 3 further options are
%               - 'xy', which returns all voxels of boundary faces with
%                       fixed z-value
%               - 'xz', which returns all voxels of boundary faces with
%                       fixed y-value
%               - 'yz', which returns all voxels of boundary faces with
%                       fixed x-value
%               - 'x',  which returns all edges in x-direction
%               - 'y',  which returns all edges in y-direction
%               - 'z',  which returns all edges in z-direction
%--------------------------------------------------------------------------
% OUTPUT
%   bndry   if 'type' is a string, bndry is a logical T_1 x ... x T_D
%           array having 1 whenever the point is part of the boundary (or 
%           the chosen subpart).
%           If 'type' is a vector of strings bndry is a structure having 
%           as fields the chosen types, which contain the bndry array.
%   weights  if 'type' is a string, weights is an T_1 x ... x T_D array
%            with weights for integration along the boundary (or the chosen
%            subpart ). A trapozoidal rule on the recangules is assumed.
%            If 'type' is a vector of strings weights instead  is a
%            structure having as fields the chosen types, which contain
%            appropriate weights array.
%--------------------------------------------------------------------------
% DEVELOPER TODOS: 
% This function doesn't work for column vectors atm, need to fix it
% E.g: 
% mask = true(D,1), bndry_voxels( mask, "full" )
%--------------------------------------------------------------------------
% EXAMPLES
% %% %% D = 1 
% %% Create a mask and show it
% Sig = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
% mask = logical( Sig > 0.02 );
% mask = mask(50,:)';
% plot( mask ),
% title( 'mask' )
% clear Sig
% 
% % Note that in 1D there is only the "full" option
% bdry = bndry_voxels( mask, "full" );
% figure(1), clf,
% plot( mask + bdry )
% title( 'mask + mask of boundary'  )
% 
% %% %% Test section D = 2
% %% 2D examples - full - example demonstrating behaviour for all ones mask
% mask = ones(10,10);
% [bndry, weights] = bndry_voxels( logical(mask), "full" )
% 
% %% 2D examples - sides
% [bndry, weights] = bndry_voxels( logical(mask), "x" )
% [bndry, weights] = bndry_voxels( logical(mask), "y" )
% 
% %% Simple box example
% % Create a simple mask and show it
% mask = zeros( [ 15 15 ] );
% mask( 4:12, 5:10 ) = 1;
% mask = logical( mask );
% figure(1), clf,
% imagesc( mask(:,:) ),
% colorbar
% title( 'mask' )
% 
% %% % plot the different options of boundary
% %% Fixed x option
% bdry = bndry_voxels( mask, "x" );
% figure(2), clf,
% imagesc( mask + bdry ), colorbar
% title("boundary for fixed x directions")
% 
% %% Fixed y option
% bdry = bndry_voxels( mask, "y" );
% figure(3), clf,
% imagesc( mask + bdry ), colorbar
% title("boundary for fixed y directions")
% 
% %% Fixed "full" option
% bdry = bndry_voxels( mask, "full" );
% figure(4), clf,
% imagesc( mask + bdry ), colorbar
% title("all boundary points")
% 
% %% Example demonstrating the bdry fix
% % Create a simple mask and show it
% mask = true( [ 10 10 ] );
% mask( 1:3, 1:4 ) = 0;
% mask = logical( mask );
% figure(1), clf,
% imagesc( mask(:,:) ),
% colorbar
% title( 'mask' )
% 
% %% % plot the different options of boundary
% %% Fixed y option
% bdry = bndry_voxels( mask, "y" );
% figure(2), clf,
% imagesc( mask + bdry ), colorbar
% title("boundary for fixed y directions")
% 
% %% Fixed x option
% bdry = bndry_voxels( mask, "x" );
% figure(3), clf,
% imagesc( mask + bdry ), colorbar
% title("boundary for fixed x directions")
% 
% %% Fixed "full" option
% bdry = bndry_voxels( mask, "full" );
% figure(4), clf,
% imagesc( mask + bdry ), colorbar
% title("all boundary points")
% 
% %% Complicated mask example
% % Generate mask
% Sig = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
% mask = logical( Sig > 0.02 & Sig < 1.1 );
% figure(1), clf,
% imagesc( mask ), colorbar,
% title("mask")
% clear Sig
% 
% %% % plot the different options of boundary
% %% Fixed y option
% bdry = bndry_voxels( mask, "y" );
% figure(2), clf,
% imagesc( mask + bdry ), colorbar
% title("boundary for fixed y directions")
% 
% %% Fixed x option
% bdry = bndry_voxels( mask, "x" );
% figure(3), clf,
% imagesc( mask + bdry ), colorbar
% title("boundary for fixed x directions")
% 
% %% Fixed "full" option
% bdry = bndry_voxels( mask, "full" );
% figure(4), clf,
% imagesc( mask + bdry ), colorbar
% title("all boundary points")
% 
% %% Fixed x+y value
% xbdry = bndry_voxels( mask, "x" );
% ybdry = bndry_voxels( mask, "y" );
% figure(3), clf,
% imagesc( (ybdry & xbdry) + xbdry ), colorbar
% title("boundary for fixed x directions")
% 
% %% %% Test section D = 3
% %% % Simple usage example on a box
% % Create a mask and show it
% mask = zeros( [ 5 5 5 ] );
% mask( 2:4, 2:4, 2:4 ) = 1;
% figure(1), clf,
% imagesc( mask(:,:,2) ), colorbar
% title("slice of mask")
% 
% % Get all boundary voxels seperate by type
% [ bdry, weights ] = bndry_voxels( logical( mask ) )
% 
% %% % Example demonstrating its use for high resolution masks
% % Create a mask and show it
% mask = true( [ 3 3 3 ] );
% mask( 1, 1, 1 ) = 0;
% visualize_bndry3D( mask, 1, ["x"], 40, [ -102 16 ] )
% title("Resolution increased mask + x-edges")
% 
% resadd = 1;
% mask_hr = mask_highres( mask, resadd, ceil(resadd/2) );
% 
% % Get all boundary voxels seperate by type
% [ bdry_xy, weights_xy ] = bndry_voxels( logical( mask_hr ) )
% 
% % Visualize all edges
% visualize_bndry3D( mask, 1, [ "x", "y", "z" ], 40, [ -102 16 ] )
% title("All edge points from bndry_voxels")
% 
% % Visualize all edges
% visualize_bndry3D( mask, 1, [ "xy", "yz", "xz" ], 40, [ -102 16 ] )
% title("All faces points from bndry_voxels")
%--------------------------------------------------------------------------
% AUTHORS: Fabian Telschow, Samuel Davenport
%--------------------------------------------------------------------------

%% Check input and get important constants from the mandatory input
%--------------------------------------------------------------------------
% Check whether the mask is logical
if ~islogical( mask )
    error( "The mask must be a logical array!" );
end

% Get the size of the mask
s_mask = size( mask );

% Get the dimension
D = length( s_mask );
if D == 2 && ( any( s_mask == 1 ) )
    D = 1;
end

% Make a larger image so that masked voxels at the boundary of the image
% will be judged to be on the boundary
[ larger_image, locs ] = pad_vals( mask );

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'type', 'var' )
   switch D
       case 2
           type = [ "x", "y" ];
       case 3
           type = [ "x", "y", "z", "xy", "xz", "yz" ];
   end
end

if ischar(type)
    type = string(type); 
    % @Fabian: Without this it seems that a structure is returned is that 
    % intended? E.g. if the input is 'full' rather than "full". To see this
    % behaviour comment this if clause out.
    % @Sam: Design idea was to be compatible with what we had before, i.e.
    % if a single type is selected always return an array. However, as soon
    % as it is a vector of types it should return a structure. As far as I
    % can see it does exactly that and that's why I needed the line.
elseif ~isstring( type )
    error( "'type' input must be a string." )
end

%% Main function
%--------------------------------------------------------------------------
% Preallocate the structures for the output
bndry   = struct();
weights = struct();
angle   = struct();
orientation = struct();

% Compute the full boundary
bndry.full = ( dilate_mask( ~larger_image, 1 ) ) & larger_image;
% bndry.full = logical( imdilate( ~larger_image, ones( ones(1, D) * 3 ) ) ) & ... larger_image;
                    
% Compute parts of the boundary if neccessary
if D == 2 || D == 3
    switch D
        case 2
            % Compute the y-boundaries
            ybndry = logical( imdilate( ~bndry.full,...
                                [ [0 0 0]; [1 1 1]; [0 0 0] ] ) ) & ...
                                bndry.full;
            ybndry = ybndry( locs{:} );
            
            % Compute the x-boundaries            
            xbndry = logical( imdilate( ~bndry.full,...
                                [ [0 1 0]; [0 1 0]; [0 1 0] ] ) ) & ...
                                    bndry.full;
            xbndry = xbndry( locs{:} );
            
            % Get full weights
            weights.full = inf2zero( 1./ ( xbndry + ybndry ) );
            
            % Check whether vertical boundary parts should be outputed
            if any( strcmp( type, "y" ) )
                bndry.y   = ybndry;
                weights.y = weights.full;
                weights.y( ~bndry.y ) = 0;
            end
            
            % Check whether horizontal boundary parts should be outputed
            if any( strcmp( type, "x" ) )
                bndry.x   = xbndry;
                weights.x = weights.full;
                weights.x( ~bndry.x ) = 0;
            end
            
            % Cut the full bndry down to the correct size
            bndry.full = bndry.full( locs{:} );
            
            if ~( any( strcmp( type, "x" ) )...
                   || any( strcmp( type, "y" ) )...
                   || any( strcmp( type, "full" ) ) )
                error( "Version must include 'full', 'x' or 'y'." );
            end
            
        case 3            
            if any( strcmp( type, "xy" ) ) || any( strcmp( type, "x" ) )...
                    || any( strcmp( type, "y" ) ) || any( strcmp( type, "z" ) )
                % Dilate kernel
                h = zeros( ones(1, D) * 3 );
                h( :, :, 2 ) = 1;
                
                % First dilation removes the edges, which do not belong to
                % a face in the specified plane. Second dilation recovers
                % the edges in the plane which belonged to a face and have
                % been removed
                bndry.xy = imdilate( ~logical( imdilate( ~bndry.full, h ) ), h );
                
                % Get logical field to obtain the orientation array
                maskp = logical( pad_vals( mask ) );
                tmp = zeros( s_mask + 2 );
                tmp( :, :, 1:end-1 ) = bndry.xy( :, :, 1:end-1) &...
                                                       maskp( :, :, 2:end);
                tmp = tmp( locs{:} );
                
                % Reduce boundary to original size
                bndry.xy = bndry.xy( locs{:} );
                
                % Orientation 1 if euclidean normal points outwards, -1 else
                orientation.xy = zeros( s_mask );
                orientation.xy( tmp == 0 ) = 1;
                orientation.xy( tmp == 1 ) = -1;
                orientation.xy( ~bndry.xy ) = 0;
                
                % Preallocate the weights array
                weights.xy = zeros( s_mask );
                % Find weights for integration
                for z = 1:s_mask(3)
                    weights.xy( :, :, z ) = getweights( bndry.xy( :, :, z ) );
                end
                weights.xy = weights.xy;
            end
                
            if any( strcmp( type, "xz" ) ) || any( strcmp( type, "x" ) )...
                    || any( strcmp( type, "y" ) ) || any( strcmp( type, "z" ) )
                % Dilate kernel
                h = zeros( ones(1, D) * 3 );
                h( :, 2, : ) = 1;
                
                % First dilation removes the edges, which do not belong to
                % a face in the specified plane. Second dilation recovers
                % the edges in the plane which belonged to a face and have
                % been removed
                bndry.xz = imdilate( ~logical( imdilate( ~bndry.full, h ) ), h );
                
                % Get logical field to obtain the orientation array
                maskp = logical( pad_vals( mask ) );
                tmp = zeros( s_mask + 2 );
                tmp( :, 1:end-1, : ) = bndry.xz( :, 1:end-1, :) &...
                                                       maskp( :, 2:end, : );
                tmp = tmp( locs{:} );
                
                % Reduce boundary to original size
                bndry.xz = bndry.xz( locs{:} );
                
                % Orientation 1 if euclidean normal points outwards, -1 else
                orientation.xz = zeros( s_mask );
                orientation.xz( tmp == 0  ) = 1;
                orientation.xz( tmp == 1  ) = -1;
                orientation.xz( ~bndry.xz ) = 0;
                
                % Preallocate the weights array
                weights.xz = zeros( s_mask );
                % Find weights for integration
                for y = 1:s_mask(2)
                    weights.xz( :, y, : ) = getweights( squeeze(...
                                                    bndry.xz( :, y, : ) ) );
                end
                weights.xz = weights.xz;             
            end
                
            if any( strcmp( type, "yz" ) ) || any( strcmp( type, "y" ) )...
                    || any( strcmp( type, "z" ) ) || any( strcmp( type, "x" ) )
                % Dilate kernel
                h = zeros( ones(1, D) * 3 );
                h( 2, :, : ) = 1;
                
                % First dilation removes the edges, which do not belong to
                % a face in the specified plane. Second dilation recovers
                % the edges in the plane which belonged to a face and have
                % been removed
                bndry.yz = imdilate( ~logical( imdilate( ~bndry.full, h ) ), h );

                % Get logical field to obtain the orientation array
                maskp = logical( pad_vals( mask ) );
                tmp = zeros( s_mask + 2 );
                tmp( 1:end-1, :, : ) = bndry.yz( 1:end-1, :, : ) &...
                                                       maskp( 2:end, :, : );
                tmp = tmp( locs{:} );
                
                % Reduce boundary to original size
                bndry.yz = bndry.yz( locs{:} );
                
                % Orientation 1 if euclidean normal points outwards, -1 else
                orientation.yz = zeros( s_mask );
                orientation.yz( tmp == 0  ) = 1;
                orientation.yz( tmp == 1  ) = -1;
                orientation.yz( ~bndry.yz ) = 0;
                
                % Preallocate the weights array
                weights.yz = zeros( s_mask );
                % Find weights for integration
                for x = 1:s_mask(1)
                    weights.yz( x, :, : ) = getweights( squeeze(...
                                                    bndry.yz( x, :, : ) ) );
                end
                weights.yz = weights.yz;
            end
            
            % Cut full weights and boundary down to correct size
            bndry.full = bndry.full( locs{:} );
            
            % Get overall weights
            weights.full = getweights( mask );
            
            if any( strcmp( type, "x" ) )
                % Get the x edges
                bndry.x = bndry.xy & bndry.xz;
 
                % Get the angle of the edge in euclidean space.
                % Note that the openeing angle is either 3/4 or 1/4
                angle.x = zeros( s_mask + 2);
                for x = 1:(s_mask(1)+2)
                    angle.x( x, :, : ) = getweights( squeeze(...
                                                    larger_image( x, :, : ) ) );
                end
                angle.x = angle.x( locs{:} );
                angle.x( ~bndry.x ) = 0;
                angle.x( angle.x == 0.5 | angle.x == 0.25 ) = 2 * pi / 4;
                angle.x( angle.x == 1   | angle.x == 0.75 ) = 2 * pi * 3 / 4;

                % Get the integration weights
                weights.x = zeros( s_mask );
                weights.x( bndry.x ) = 1;
                % Find the corners of each edge and halft its weight
                weights.x( bndry.xy & bndry.xz & bndry.yz ) = ...
                weights.x( bndry.xy & bndry.xz & bndry.yz ) / 2;
            end
            
            if any( strcmp( type, "y" ) )
                % Get the y edges
                bndry.y = bndry.xy & bndry.yz;
                
                % Get the angle of the edge in euclidean space.
                % Note that the openeing angle is either 3/4 or 1/4
                angle.y = zeros( s_mask + 2);
                for y = 1:( s_mask(1) + 2 )
                    angle.y( :, y, : ) = getweights( squeeze(...
                                                    larger_image( :, y, : ) ) );
                end
                angle.y = angle.y( locs{:} );
                angle.y( ~bndry.y ) = 0;
                angle.y( angle.y == 0.5 | angle.y == 0.25 ) = 2 * pi / 4;
                angle.y( angle.y == 1   | angle.y == 0.75 ) = 2 * pi * 3 / 4;
                
                % Get the integration weights
                weights.y = zeros( s_mask );
                weights.y( bndry.y ) = 1;              
                % Find the corners of each edge and halft its weight
                weights.y( bndry.xy & bndry.xz & bndry.yz ) = ...
                weights.y( bndry.xy & bndry.xz & bndry.yz ) / 2;  
            end
            
            if any( strcmp( type, "z" ) )
                % Get the z edges
                bndry.z = bndry.yz & bndry.xz;
                 
                % Get the angle of the edge in euclidean space.
                % Note that the openeing angle is either 3/4 or 1/4
                angle.z = zeros( s_mask + 2);
                for z = 1:(s_mask(1)+2)
                    angle.z( :, :, z ) = getweights( squeeze(...
                                                larger_image( :, :, z ) ) );
                end
                angle.z = angle.z( locs{:} );
                angle.z( ~bndry.z ) = 0;
                angle.z( angle.z == 0.5 | angle.z == 0.25 ) = 2 * pi / 4;
                angle.z( angle.z == 1   | angle.z == 0.75 ) = 2 * pi * 3 / 4;

                % Get the integration weights
                weights.z = zeros( s_mask );
                weights.z( bndry.z ) = 1;
                % Find the corners of each edge and halft its weight
                weights.z( bndry.xy & bndry.xz & bndry.yz ) = ...
                weights.z( bndry.xy & bndry.xz & bndry.yz ) / 2;
            end
            
            % Check whether the user entered at least on appropriate
            % type input
            if ~( any( strcmp( type, "full" ) )...
                    || any( strcmp( type, "xy" ) )...
                    || any( strcmp( type, "xz" ) )...
                    || any( strcmp( type, "yz" ) )...
                    || any( strcmp( type, "x" ) )...
                    || any( strcmp( type, "y" ) )...
                    || any( strcmp( type, "z" ) ) )
                error( "'Type' must be chosen according to description." );
            end     
    end
else
    bndry.full = bndry.full( locs{:} );
    if D > 3
        weights.full = getweights( mask );
    else
        weights.full = double( mask );
        weights.full( bndry.full ) = 1/2;
    end
end

% Make output arrays if type is a string
if length( type ) == 1
    bndry    = bndry.(type);
    weights  = weights.(type);
elseif ~any( strcmp( type, "full" ) )
    bndry = rmfield( bndry, 'full' );
    weights = rmfield( weights, 'full' );
end

return