function cov_array = vectcov( array1, array2, dimension, normalize )
% VECTCOV( array1, array2, dimension, normalize ) calculates the covariance
% of two arrays along a given dimension.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  array1      a D dimensional array of data
%  array2      a D dimensional array of data with the same size as array1
% Optional 
%  dimension   the dimension along which to compute the covariance.
%              Default, is to take dimension to be D
%  normalize   0/1 whether or not to normalize the covariance by N/(N-1)
%--------------------------------------------------------------------------
% OUTPUT
%  cov_array   a D-1 dimensional array with the same size as array1 and
%              array2 except for the specified dimension where each entry 
%              isthe covariance of points computed along the specified
%              dimension
%--------------------------------------------------------------------------
% EXAMPLES
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Obtain the array sizes 
s1 = size(array1); s2 = size(array2);

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist('dimension', 'var')
    % default option of dimension
    dimension = length(s1);
end
if ~exist('normalize', 'var')
    % default option of normalize
    normalize = 1;
end

%% Error checking
%--------------------------------------------------------------------------
if ~isequal(s1, s2)
   error('The arrays must have the same dimension')
end
if dimension > length(s1)
    error('The dimension must be <= the number of dimensions of the arrays')
end

%% Main function
%--------------------------------------------------------------------------
% Calculate the mean of the first araay along the specified dimension
array1_mean = mean(array1, dimension);

% Calculate the covariance between array1 and array2 along the specified
% dimension
cov_array = mean((array1 - array1_mean).*array2, dimension);

% Compute the total number of subjects
N = size(array1, dimension);

% If the normalize option is specified scale the covariance by the
% normalizing factor
if normalize == 1
    cov_array = (cov_array)*(N/(N-1));
end

end
