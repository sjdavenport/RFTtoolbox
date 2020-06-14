function out = convind( ind, size_of_array, conv2what)
% CONVIND( ind, size_of_array, conv2what ) switches between vector indexing and 
% coordinate indexing. 
% For a size vector x = size [x_1, ..., x_n] and an array A of size x we can 
% convert A to a vector of length x_1*...*x_n via A(:). Similarly for any
% vector v of length x_1*...*x_n we can convert it to an array of size x via
% reshape(v, x). The location of an element in the array is given by a
% vector of length n, and this corresponds to a unique number between 1 and  
% x_1*...*x_n in the vector representation. This function converts between
% the two.
%--------------------------------------------------------------------------
% ARGUMENTS
% ind       A numeric vector. If its of length one it corresponds to the
%           location in the vector representation and if its of length 
%           greater than one it corresponds to the coordinate vector.
% conv2what       0/1/2, converts the indices to fsleyes indices if 1 and
%                 to aligned anatomical coordinates (for use in papers) if 
%                 set to 2. Just returns the coordinates (the default) is
%                 affected if set to 0. 
% size_of_array   A numeric vector. This corresponds to [x_1, ..., x_n] in 
%                 the above description. This is defaulted to [91,109,91]
%                 which is the default size of the 2mm MNI brain image 
%                 in voxels.
%--------------------------------------------------------------------------
% OUTPUT
% out       out is a coordinate vector if ind has length 1 and it is a
%           number if ind is a vector of length greater than 1.
%--------------------------------------------------------------------------
% EXAMPLES
% convind([32,15,37])
% convind(358390)
%
% convind([13,15], [250,250])
% convind(3513, [250,250])
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if nargin < 2
    stdsize = [91,109,91];
    size_of_array = stdsize;
end
if nargin < 3
    conv2what = 0;
end

len_ind = length(ind);
% warning('if there is an error oculd be because you switched the 2nd two arguments')

D = length(size_of_array);
if len_ind > 1 && ~(len_ind == D)
   error('The length of ind must be 1 (vector representation) or be the same length as the size of the array (array representation).') 
end

if D == 2 && (size_of_array(1) == 1 || size_of_array(1) == 1)
    out = ind; %I.e. as in 1D there's no need to change anything!
    return
end

if len_ind == 1
    %Vector to coordinate change
    if D == 3
        [x,y,z] = ind2sub(size_of_array, ind);
        out = [x,y,z];
    elseif D == 2
        [x,y] = ind2sub(size_of_array, ind);
        out = [x,y];
    end
    if conv2what == 1
        out = out - 1;
    elseif conv2what == 2
        out = out - 1;
        out(1) = 90 - 2*out(1);
        out(2) = 2*out(2) - 126;
        out(3) = 2*out(3) - 72;
    end 
else
    %Coordinate to vector change
    ind_string = sprintf('%.0f,', ind);
    ind_string = ind_string(1:end-1);
    out = eval(strcat('sub2ind(size_of_array, ',ind_string,')'));
end

out = squeeze(out);

end

