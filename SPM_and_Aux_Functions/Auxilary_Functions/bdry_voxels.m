function bdry_voxels = bdry_voxels( mask, version )
% This function computes a high resolution version of a given mask.
% It has the option to enlarge the mask region by resAdd to use shifted
% boundaries in LKC estimation. This is required in the interpretation of
% values at voxels as the center values of rectangular domains. 
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   mask    an logical T_1 x ... x T_D array.
%   resAdd  the amount of equidistant voxels introduced inbetween the
%           voxels  
%--------------------------------------------------------------------------
% OUTPUT
%   mask_hr     an logical hr_T_1 x ... x hr_T_D array. Here hr_T_i = 
%   plot_switch logical to show educational plots explaining the code of
%               the function. Default is 0, i.e. no plots.
%--------------------------------------------------------------------------
% DEVELOPER TODOs:
% -------------------------------------------------------------------------
% AUTHORS: Fabian Telschow
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% main function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if version == "full"
    bdry_voxels = logical( imdilate( ~mask, ones( ones(1, D)*3 ) ) ) & ...
                            mask;
elseif D==2 || D==3

    switch D
        case 2
            if version == "y"
                bdry_voxels = logical( imdilate( ~mask,...
                                    [[0 0 0];[1 1 1];[0 0 0]] ) ) & ...
                                        mask;
            elseif version == "x"
                bdry_voxels = logical( imdilate( ~mask,...
                                    [[0 1 0];[0 1 0];[0 1 0]] ) ) & ...
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
                bdry_voxels = logical( imdilate( ~mask, h ) ) &  mask;
            elseif version == "yz"
                h( 2, 1, 2 ) = 1;
                h( 2, 3, 2 ) = 1;
                bdry_voxels = logical( imdilate( ~mask, h ) ) &  mask;
            elseif version == "xz"
                h( 1,2,  2 ) = 1;
                h( 3,2,  2 ) = 1;
                bdry_voxels = logical( imdilate( ~mask, h ) ) &  mask;
            else
                error( "Version must be either 'full', 'x' or 'y'" );
            end
            
    end
else
    error( strcat( "For D not equal to 2 or 3 only the full",...
                   "boundary estimate is implemented" ) );
end


return