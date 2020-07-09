function [ out, array_diff ] = same_array( array1, array2, tol )
% SAME_ARRAY( array1, array2, tol )
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  array1       an array of data
%  array2       another array of data with the same dimensions as array1
% Optional
%  tol 
%--------------------------------------------------------------------------
% OUTPUT
% out           0/1 whether the arrays are the same
%--------------------------------------------------------------------------
% EXAMPLES
% array_1 = randn([3,4]); array_2 = randn([3,4]);
% same_array(array_1, array_2)
% same_array(array_1, array_1)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'tol', 'var' )
   % default option of opt1
   tol = 10^(-7);
end

%%  Main Function Loop
%--------------------------------------------------------------------------
notnan1 = ~isnan(array1);
notnan2 = ~isnan(array2);

array_diff = sum(abs(array1(notnan1) - array2(notnan2)));
if array_diff < tol
    out = 1;
else
    out = 0;
end

end

