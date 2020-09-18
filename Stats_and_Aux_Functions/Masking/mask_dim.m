function D = mask_dim( mask )
% MASK_DIM( mask ) obtains the number of dimensions of the mask
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% mask      an array of data
%--------------------------------------------------------------------------
% OUTPUT
% D         the number of dimensions of the mask
%--------------------------------------------------------------------------
% EXAMPLES
% D = mask_dim(ones(10,1))
% D = mask_dim(ones(100,100))
% D = mask_dim(ones(100,100,100))
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

% Find the size of the mask
s_mask = size(mask);

% Find the length of the size vector
D = length(s_mask);

if D == 2
    % If the mask is of the form [nvox,1] or [1,nvox] then D = 1!
   if s_mask(1) == 1 ||  s_mask(2) == 1
       D = 1;
   end
end

end
