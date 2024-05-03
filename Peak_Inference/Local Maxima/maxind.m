function max_index = maxind( map, ascoord )
% MAX_IND( map ) find the index corresponding to the maximum of the 
% 91by109by91 map.
%--------------------------------------------------------------------------
% ARGUMENTS
% map       a 91by109by91 array.
% ascoord   0/1. Default 1. If 1 then the maximum coorindate is returned in
%           coordinate form. If 0 its represented in vector coordinate
%           form.
%--------------------------------------------------------------------------
% OUTPUT
% index     a length 3 vector with the coordinates of the maximum.
%--------------------------------------------------------------------------
% EXAMPLES
% a = zeros([91,109,91]);
% a(16,100,40) = 5;
% max_index = max_ind(a) %need to use 'all' here as (16,100,40)
% %doesn't lie within the MNI mask of the brain.
% a(max_index(1), max_index(2), max_index(3));
%--------------------------------------------------------------------------
% AUTHOR: Sam Davenport
if nargin < 2
    ascoord = 1;
end

max_val = max(map(:));
max_loc = find(map(:) == max_val);

if ascoord
    max_index = convind(max_loc);
end

end

