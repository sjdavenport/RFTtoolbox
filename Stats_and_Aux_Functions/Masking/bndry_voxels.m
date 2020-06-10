function [ bdry, weights ] = bndry_voxels( mask, version )
% This function computes a high resolution version of a given mask.
% It has the option to enlarge the mask region by resAdd to use shifted
% boundaries in LKC estimation. This is required in the interpretation of
% values at voxels as the center values of rectangular domains. 
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   mask     a logical T_1 x ... x T_D array.
%   version  a string indicating which part of the boundary is obtained.
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
%   bdry     a logical T_1 x ... x T_D array having 1 whenever the point
%            is part of the boundary (or the chosen subpart).
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
% bdry = bdry_voxels( mask, "full" );
% figure(1), clf,
% plot( mask + bdry )
% title( 'mask + mask of boundary'  )
% 
% %% Test section D = 2
% % Example demonstrating behaviour for all ones mask
% % 2D examples - full
% mask = ones(10,10);
% bdry_voxels( logical(mask), 'full' )
% 
% % 2D examples - sides
% bdry_voxels( logical(mask), 'x' )
% bdry_voxels( logical(mask), 'y' )
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
% bdry = bdry_voxels( mask, "x" );
% figure(2), clf,
% imagesc( mask + bdry ), colorbar
% title("boundary for fixed y directions")
% 
% % fixed x option
% bdry = bdry_voxels( mask, "y" );
% figure(3), clf,
% imagesc( mask + bdry ), colorbar
% title("boundary for fixed x directions")
% 
% % fixed "full" option
% bdry = bdry_voxels( mask, "full" );
% figure(4), clf,
% imagesc( mask + bdry ), colorbar
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
% bdry = bdry_voxels( mask, "x" );
% figure(2), clf,
% imagesc( mask + bdry ), colorbar
% title("boundary for fixed y directions")
% 
% % fixed x option
% bdry = bdry_voxels( mask, "y" );
% figure(3), clf,
% imagesc( mask + bdry ), colorbar
% title("boundary for fixed x directions")
% 
% % fixed "full" option
% bdry = bdry_voxels( mask, "full" );
% figure(4), clf,
% imagesc( mask + bdry ), colorbar
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
% bdry = bdry_voxels( logical( mask ), "xy" )
% 
% % same as before for y
% bdry = bdry_voxels( logical( mask ), "xz" )
% 
% % same as before for x
% bdry = bdry_voxels( logical( mask ), "yz" )
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
if ~exist( 'version', 'var' )
   % default option of version
   version = 'full';
end

%% Main function
%--------------------------------------------------------------------------

% compute the full boundary
bdry_full = logical( imdilate( ~larger_image, ones( ones(1, D) * 3 ) ) ) & ...
                        larger_image;
                    
% compute parts of the boundary if neccessary
if D==2 || D==3
    switch D
        case 2
            weights = NaN * bdry_full;
            bdry = bdry_full;
            if version == "y"
                bdry = logical( imdilate( ~bdry,...
                                    [ [0 0 0]; [1 1 1]; [0 0 0] ] ) ) & ...
                                        bdry;
            elseif version == "x"
                bdry = logical( imdilate( ~bdry,...
                                    [ [0 1 0]; [0 1 0]; [0 1 0] ] ) ) & ...
                                        bdry;
            elseif version ~= "full"
                error( "Version must be either 'full', 'x' or 'y'" );
            end
            
        case 3
            h = zeros( ones(1, D) * 3 );
            % First dilation removes the edges, which do not belong to a
            % face in the specified plane. Second dilation recovers the
            % edges in the plane which belonged to a face and have been
            % removed
            if version == "xy"
                h( :, :, 2 ) = 1;
                bdry = imdilate( ~logical( imdilate( ~bdry_full, h ) ), h );
                
                % modify filter to count orthogonal neighbours, since they
                % define to how many faces the bdry voxel belongs to
                h( :, :, 2) = 0;
                h( 2, :, 2) = 1;
                h( :, 2, 2) = 1;
                h( 2, 2, 2) = 0;
                % count neighbours in the surface
                weights =  convn( bdry, h, 'same' );
                
                % set the weights correctly
                weights( ~bdry ) = 0;
                weights( weights == 2 ) = 1/4;
                weights( weights == 3 ) = 2/4;
                weights( weights == 4 ) = 1;
                
            elseif version == "yz"
                h( :, 2, : ) = 1;
                bdry = imdilate( ~logical( imdilate( ~bdry_full, h ) ), h );
                
                % modify filter to count orthogonal neighbours, since they
                % define to how many faces the bdry voxel belongs to
                h( :, 2, :) = 0;
                h( 2, 2, :) = 1;
                h( :, 2, 2) = 1;
                h( 2, 2, 2) = 0;
                % count neighbours in the surface
                weights =  convn( bdry, h, 'same' );
                
                % set the weights correctly
                weights( ~bdry ) = 0;
                weights( weights == 2 ) = 1/4;
                weights( weights == 3 ) = 2/4;
                weights( weights == 4 ) = 1;
                
            elseif version == "xz"
                h( 2, :, : ) = 1;
                bdry = imdilate( ~logical( imdilate( ~bdry_full, h ) ), h );
                
                % modify filter to count orthogonal neighbours, since they
                % define to how many faces the bdry voxel belongs to
                h( 2, :, :) = 0;
                h( 2, 2, :) = 1;
                h( 2, :, 2) = 1;
                h( 2, 2, 2) = 0;
                % count neighbours in the surface
                weights =  convn( bdry, h, 'same' );
                
                % set the weights correctly
                weights( ~bdry ) = 0;
                weights( weights == 2 ) = 1/4;
                weights( weights == 3 ) = 2/4;
                weights( weights == 4 ) = 1;
                
            elseif version == "full"
                bdry = bdry_full;
                weights = NaN * bdry;
            else
                error( "Version must be either 'full', 'x' or 'y'" );
            end
            
    end
else
    error( strcat( "For D not equal to 2 or 3 only the full",...
                   "boundary estimate is implemented" ) );
end

% Remove the outer voxels
bdry    = bdry( locs{:} );
weights = weights( locs{:} );

return