function [ Dim, D ] = getdim( mask )
% GETDIM( mask ) obtains the dimensions from the mask including in 1D when
% the mask is of the form [nvox,1] or [1,nvox].
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  mask       a 0/1 mask of the data
%--------------------------------------------------------------------------
% OUTPUT
%  Dim   the dimensions of the mask (note there is a special
%        condition allowing dimensional 1 data to be set)
%  D     the number of dimensions
%--------------------------------------------------------------------------
% EXAMPLES
% %% 1D
% [ Dim, D ] = getdim(ones(10,1))
% 
% %% 2D
% [ Dim, D ] = getdim(ones(10,10))
% 
% %% 3D
% [ Dim, D ] = getdim(ones(10,10,10))
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
s_mask = size(mask);

%%  Main Function Loop
%--------------------------------------------------------------------------
if length(s_mask) == 2 && (s_mask(1) == 1 || s_mask(2) == 1)
    % Dimension 1
    Dim = s_mask(s_mask > 1);
    D = 1;
else
    % Other dimensions
    Dim = size(mask);
    D = length(Dim);
end

end

