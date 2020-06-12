function [ bndry, weights ] = bndry_voxels( mask, type )
% This function computes a high resolution version of a given mask.
% It has the option to enlarge the mask region by resAdd to use shifted
% boundaries in LKC estimation. This is required in the interpretation of
% values at voxels as the center values of rectangular domains. 
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
%--------------------------------------------------------------------------
% OUTPUT
%   bndry    if 'type' is a string, bndry is a logical T_1 x ... x T_D
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
% EXAMPLES
% %% Test section D = 1
% % create a mask and show it
% Sig = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
% mask = logical( Sig > 0.02 );
% mask = mask(50,:)';
% plot( mask ),
% title( 'mask' )
% clear Sig
% 
% % note in 1D there is only the 'full' option
% bndry = bndry_voxels( mask, "full" );
% figure(1), clf,
% plot( mask + bndry )
% title( 'mask + mask of boundary'  )
% 
% %% Test section D = 2
% % Example demonstrating behaviour for all ones mask
% % 2D examples - full
% mask = ones(10,10);
% bndry_voxels( logical(mask), 'full' )
% 
% % 2D examples - sides
% bndry_voxels( logical(mask), 'x' )
% bndry_voxels( logical(mask), 'y' )
% 
% % simple box example
% close all
% % create a simple mask and show it
% mask = zeros( [ 15 15 ] );
% mask( 4:12, 5:10 ) = 1;
% mask = logical( mask );
% figure(1), clf,
% imagesc( mask(:,:) ),
% colorbar
% title( 'mask' )
% 
% %%% plot the different options of boundary
% % fixed y option
% bndry = bndry_voxels( mask, "x" );
% figure(2), clf,
% imagesc( mask + bndry ), colorbar
% title("boundary for fixed y directions")
% 
% % fixed x option
% bndry = bndry_voxels( mask, "y" );
% figure(3), clf,
% imagesc( mask + bndry ), colorbar
% title("boundary for fixed x directions")
% 
% % fixed "full" option
% bndry = bndry_voxels( mask, "full" );
% figure(4), clf,
% imagesc( mask + bndry ), colorbar
% title("all boundary points")
% 
% %% complicated mask example
% % generate mask
% Sig = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
% mask = logical( Sig > 0.02 & Sig < 1.1 );
% figure(1), clf,
% imagesc( mask ), colorbar,
% title("mask")
% clear Sig
% 
% %% plot the different options of boundary
% % fixed y option
% bndry = bndry_voxels( mask, "x" );
% figure(2), clf,
% imagesc( mask + bndry ), colorbar
% title("boundary for fixed y directions")
% 
% % fixed x option
% bndry = bndry_voxels( mask, "y" );
% figure(3), clf,
% imagesc( mask + bndry ), colorbar
% title("boundary for fixed x directions")
% 
% % fixed "full" option
% bndry = bndry_voxels( mask, "full" );
% figure(4), clf,
% imagesc( mask + bndry ), colorbar
% title("all boundary points")
% 
% %% Test section D = 3
% % create a mask and show it
% mask = zeros( [ 5 5 5 ] );
% mask( 2:4, 2:4, 2:4 ) = 1;
% figure(1), clf,
% imagesc( mask(:,:,2) ), colorbar
% title("slice of mask")
% 
% % get boundary voxels lying in z value planes or better having a not to the
% % mask connected face pointing into the z direction
% bndry = bndry_voxels( logical( mask ), "xy" )
% 
% % same as before for y
% bndry = bndry_voxels( logical( mask ), "xz" )
% 
% % same as before for x
% bndry = bndry_voxels( logical( mask ), "yz" )
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
if D == 2 && s_mask(2) == 1
    D = 1;
end

% Make a larger image so that masked voxels at the boundary of the image
% will be judged to be on the boundary
[ larger_image, locs ] = pad_vals( mask );

%%  add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'type', 'var' )
   switch D
       case 2
           type = ["x", "y"];
       case 3
           type = ["x", "y", "z", "xy", "xz", "yz"];
   end
end

if ~isstring( type )
    error( "'type' input must be a string." )
end

%% Main function
%--------------------------------------------------------------------------
% Preallocate the structures for the output
bndry   = struct();
weights = struct();

% compute the full boundary
bndry.full = logical( imdilate( ~larger_image, ones( ones(1, D) * 3 ) ) ) & ...
                        larger_image;
                    
% compute parts of the boundary if neccessary
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
            
            % check whether vertical boundary parts should be outputed
            if any( strcmp( type, "y" ) )
                bndry.y   = ybndry;
                weights.y = weights.full;
                weights.y( ~bndry.y ) = 0;
            end
            
            % check whether horizontal boundary parts should be outputed
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
                    || any( strcmp( type, "y" ) ) || any( strcmp( type, "x" ) )
                % Dilate kernel
                h = zeros( ones(1, D) * 3 );
                h( :, :, 2 ) = 1;
                
                % First dilation removes the edges, which do not belong to
                % a face in the specified plane. Second dilation recovers
                % the edges in the plane which belonged to a face and have
                % been removed
                bndry.xy = imdilate( ~logical( imdilate( ~bndry.full, h ) ), h );
                bndry.xy = bndry.xy( locs{:} );
                
                % Preallocate the weights array
                weights.xy = zeros( s_mask );
                % Find weights for integration
                for z = 1:s_mask(3)
                    weights.xy( :, :, z ) = getweights( bndry.xy( :, :, z ) );
                end
                weights.xy = weights.xy;
            end
                
            if any( strcmp( type, "xz" ) ) || any( strcmp( type, "x" ) )...
                    || any( strcmp( type, "z" ) ) || any( strcmp( type, "x" ) )
                % Dilate kernel
                h = zeros( ones(1, D) * 3 );
                h( :, 2, : ) = 1;
                
                % First dilation removes the edges, which do not belong to
                % a face in the specified plane. Second dilation recovers
                % the edges in the plane which belonged to a face and have
                % been removed
                bndry.xz = imdilate( ~logical( imdilate( ~bndry.full, h ) ), h );
                bndry.xz = bndry.xz( locs{:} );
                
                % Preallocate the weights array
                weights.xz = zeros( s_mask );
                % Find weights for integration
                for y = 2:s_mask(2)
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
                bndry.yz = bndry.yz( locs{:} );
                
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
 
                % Get the neighbours. Note that the openeing angle on an
                % edge is either 3/4 or 1/4
                weights.x = zeros( s_mask + 2);
                for x = 1:(s_mask(1)+2)
                    weights.x( x, :, : ) = getweights( squeeze(...
                                                    larger_image( x, :, : ) ) );
                end
                weights.x = weights.x( locs{:} );
                weights.x( ~bndry.x ) = 0;
                weights.x( weights.x == 0.5 ) = 1/4;
                weights.x( weights.x == 1   ) = 3/4;
                                 
                % Find the corners of each edge and halft its weight
                weights.x( bndry.xy & bndry.xz & bndry.yz ) = ...
                weights.x( bndry.xy & bndry.xz & bndry.yz ) / 2; 
            end
            
            if any( strcmp( type, "y" ) )
                % Get the y edges
                bndry.y = bndry.xy & bndry.yz;
                
                % Get the neighbours. Note that the openeing angle on an
                % edge is either 3/4 or 1/4
                weights.y = zeros( s_mask + 2);
                for y = 1:(s_mask(1)+2)
                    weights.y( :, y, : ) = getweights( squeeze(...
                                                    larger_image( :, y, : ) ) );
                end
                weights.y = weights.y( locs{:} );
                weights.y( ~bndry.y ) = 0;
                weights.y( weights.y == 0.5 ) = 1/4;
                weights.y( weights.y == 1   ) = 3/4;
                 
                % Find the corners of each edge and halft its weight
                weights.y( bndry.xy & bndry.xz & bndry.yz ) = ...
                weights.y( bndry.xy & bndry.xz & bndry.yz ) / 2;  
            end
            
            if any( strcmp( type, "z" ) )
                % Get the z edges
                bndry.z = bndry.yz & bndry.xz;
                 
                % Get the neighbours. Note that the openeing angle on an
                % edge is either 3/4 or 1/4
                weights.z = zeros( s_mask + 2);
                for z = 1:(s_mask(1)+2)
                    weights.z( :, :, z ) = getweights( squeeze(...
                                                larger_image( :, :, z ) ) );
                end
                weights.z = weights.z( locs{:} );
                weights.z( ~bndry.z ) = 0;
                weights.z( weights.z == 0.5 ) = 1/4;
                weights.z( weights.z == 1   ) = 3/4;
                 
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
    weights.full = weights.full( locs{:} );  
end

% Make output arrays if type is a string
if length( type ) == 1
    bndry    = bndry.(type);
    weights = weights.(type);
elseif ~any( strcmp( type, "full" ) )
    bndry = rmfield( bndry, 'full' );
    weights = rmfield( weights, 'full' );
end

return