function image = reshape3D( vector )
% RESHAPE3D( array ) converts a 91x109x91 vector to a 91by109by91 image.
%--------------------------------------------------------------------------
% ARGUMENTS
% vector    a 91x109x91 vector
%--------------------------------------------------------------------------
% OUTPUT
% image     an 91by109by91 image
%--------------------------------------------------------------------------
% AUTHOR: Sam Davenport.

if length(vector) ~= prod([91,109,91])
    error('The vector must be of length 91x109x91')
end
    
image = reshape(vector, [91,109,91]);

end

