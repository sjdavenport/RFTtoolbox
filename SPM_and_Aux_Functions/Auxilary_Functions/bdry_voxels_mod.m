function bdry = bdry_voxels( mask, version )
% This function computes a high resolution version of a given mask.
% It has the option to enlarge the mask region by resAdd to use shifted
% boundaries in LKC estimation. This is required in the interpretation of
% values at voxels as the center values of rectangular domains. 
%
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
%   bdry     an logical T_1 x ... x T_D array having 1 whenever the point
%            is part of the boundary (or the chosen subpart).
%--------------------------------------------------------------------------
% DEVELOPER TODOs:
% -------------------------------------------------------------------------
% EXAMPLES
%--------------------------------------------------------------------------
% AUTHORS: Fabian Telschow
%--------------------------------------------------------------------------

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Check input and get important constants from the mandatory input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check whether the mask is logical
if ~islogical( mask )
    error( "The mask must be a logical array!" );
end

% get the size of the mask
s_mask = size( mask );

% get the dimension
D = length( s_mask );
if D == 2 && s_mask(2) == 1
    D = 1;
end

% Make a larger image so that masked voxels at the boundary of the image
% will be judged to be on the boundary
larger_image = zeros(s_mask+2);

% Get the locations to place the inner (original) data
b = cell(1,D);
for d = 1:D
   b{d} = 2:s_mask(d)+1;
end

% Set the inner locations to be the mask
larger_image(b{:}) = mask;
mask = larger_image; %Can remove this line and rep mask with larger image

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% main function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if version == "full"
    bdry = logical( imdilate( ~mask, ones( ones(1, D) * 3 ) ) ) & ...
                            mask;
elseif D==2 || D==3

    switch D
        case 2
            if version == "y"
                bdry = logical( imdilate( ~mask,...
                                    [ [0 0 0]; [1 1 1]; [0 0 0] ] ) ) & ...
                                        mask;
            elseif version == "x"
                bdry = logical( imdilate( ~mask,...
                                    [ [0 1 0]; [0 1 0]; [0 1 0] ] ) ) & ...
                                        mask;
            else
                error( "Version must be either 'full', 'x' or 'y'" );
            end
            
        case 3
            h = zeros( ones(1, D)*3 );
            h( 2, 2, 2 ) = 1;
            if version == "xy"
                h( 2, 2, 1 ) = 1;
                h( 2, 2, 3 ) = 1;
                bdry = logical( imdilate( ~mask, h ) ) &  mask;
            elseif version == "yz"
                h( 2, 1, 2 ) = 1;
                h( 2, 3, 2 ) = 1;
                bdry = logical( imdilate( ~mask, h ) ) &  mask;
            elseif version == "xz"
                h( 1,2,  2 ) = 1;
                h( 3,2,  2 ) = 1;
                bdry = logical( imdilate( ~mask, h ) ) &  mask;
            else
                error( "Version must be either 'full', 'x' or 'y'" );
            end
            
    end
else
    error( strcat( "For D not equal to 2 or 3 only the full",...
                   "boundary estimate is implemented" ) );
end

% new_bdry = zeros(s_mask);
actual_bdry = bdry(b{:});
bdry = actual_bdry;

return