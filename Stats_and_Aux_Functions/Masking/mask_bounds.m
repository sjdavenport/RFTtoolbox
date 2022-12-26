function [ bounds, bounded_mask ] = mask_bounds( mask )
% MASK_BOUNDS finds the bounds of a mask and returns a mask without all
% the padding.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% mask - a D-dimensional array of zeros and ones
%--------------------------------------------------------------------------
% OUTPUT
% bounds - a D-by-1 cell array, where the dth entry gives the bounds of the 
%           mask in the dth direction
% bounded_mask - a logical array giving the mask restricted without the
%                extra padding
%--------------------------------------------------------------------------
% EXAMPLES
% mask = imgload('MNImask'); [ bounds, bounded_mask ] = mask_bounds( mask )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Obtain the size of the mask
s_mask = size(mask);

% Obtain the number of dimensions of the mask
D = length(s_mask);

%%  Main Function Loop
%--------------------------------------------------------------------------
bounds = cell(1, D);

for d = 1:D
    index = repmat( {':'}, 1, D );
    lower_bound_d = 0;
    slice = 0;
    
    % Obtain the lower bound in the dth direction
    while lower_bound_d == 0
        slice = slice + 1;
        index{d} = slice;
        masked_slice = mask(index{:});
        if sum(masked_slice(:)) > 0
            lower_bound_d = slice;
        end
    end
   
    %Initialize the upper_bound and the slice
    upper_bound_d = 0;
    slice = s_mask(d);
    
    % Obtain the upper bound in the dth direction
    while upper_bound_d == 0
        index{d} = slice;
        masked_slice = mask(index{:});
        if sum(masked_slice(:)) > 0
            upper_bound_d = slice;
        end
        
        % Decrease slice by 1
        slice = slice - 1;
    end
    
    bounds{d} = lower_bound_d:upper_bound_d;
end

bounded_mask = logical(mask(bounds{:}));

end

