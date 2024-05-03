function [ A, detg ] = invsym( S )
% invg( S ) calculates the determinant and the inverse of a symmetric matrix.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  S      a D dimensional array of data
%--------------------------------------------------------------------------
% OUTPUT
%  cov_array   a D-1 dimensional array with the same size as array1 and
%              array2 except for the specified dimension where each entry 
%              isthe covariance of points computed along the specified
%              dimension
%--------------------------------------------------------------------------
% EXAMPLES
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow
%--------------------------------------------------------------------------

%% Main function
%--------------------------------------------------------------------------

index  = repmat( {':'}, 1, length(size(S))-2 );

S11 = S( index{:}, 1, 1 );
S12 = S( index{:}, 1, 2 );
S13 = S( index{:}, 1, 3 );
S22 = S( index{:}, 2, 2 );
S23 = S( index{:}, 2, 3 );
S33 = S( index{:}, 3, 3 );

% Get the determinant
detg =   S11 .* S22 .* S33 + 2 * S12 .* S23 .* S13...
       - S11 .* S23.^2 - S22 .* S13.^2 - S33 .* S12.^2;

% Get the inverse by using g^-1 = adj(g)/det(g)
A   = zeros( size(S) );
A( index{:}, 1, 1 ) = S22 .* S33 - S23.^2;
A( index{:}, 1, 2 ) = S13 .* S23 - S12 .* S33;
A( index{:}, 1, 3 ) = S12 .* S23 - S13 .* S22;

A( index{:}, 2, 1 ) = A( index{:}, 1, 2 );
A( index{:}, 2, 2 ) = S11 .* S33 - S13.^2;
A( index{:}, 2, 3 ) = S12 .* S13 - S11 .* S23;

A( index{:}, 3, 1 ) = A( index{:}, 1, 3 );
A( index{:}, 3, 2 ) = A( index{:}, 2, 3 );
A( index{:}, 3, 3 ) = S22 .* S11 - S12.^2;

A = A ./ detg;

end
