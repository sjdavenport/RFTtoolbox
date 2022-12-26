function out = vech( mateorvec )
% VECH converts a matrix to a vector containing its upper triangular
% elements and converts a vector to a symmetric matrix with the upper
% triangular elements in the right place.
%--------------------------------------------------------------------------
% INPUT
%   mateorvec:  a matrix or a vector
%--------------------------------------------------------------------------
% OUTPUT
%   out:        a vector containing the upper triangular elements of the 
%               input matrix if the input is a matrix, or a symmetric matrix
%               with the upper triangular elements of the input vector in the
%               right place if the input is a vector.
%--------------------------------------------------------------------------
% EXAMPLES
%   vech([1,2;3,4])
%   vech([1,2,4])
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
smate = size(mateorvec);
if length(smate) > 2
    error('Need to enter a matrix or a vector')
end

% If a vector ensure that it is a row vector
if smate(2) == 1
    smate = flipud(smate);
end

%%  Main Function Loop
%--------------------------------------------------------------------------
if smate(1) == 1
    % (The inverse of D(D+1)/2 = smate(1));
    D = (sqrt(1+8*smate(2))-1)/2;
    
    % Initialize the matrix to return
    out = zeros(D,D);
    
    % Keep track of how many elements of the vector have been used
    len_of_vec_used = 0;
    
    % Calculate the upper triangular matrix
    for d = 1:D
        out(1:d, d) = mateorvec( (len_of_vec_used +1):(len_of_vec_used+d));
        len_of_vec_used = len_of_vec_used + d;
    end
    
    % Make the output matrix symmetric  
    out = out + tril(out',-1);
else
    % Calculate the upper triangular matrix
    upper_triangular_mate = triu(mateorvec);
    
    % Calculate the upper triangular matrix of ones
    upper_triangular_ones = triu(ones(smate));
    
    % Obtain the upper triangular matrix as a vector
    out = upper_triangular_mate(logical(upper_triangular_ones));
end

end

