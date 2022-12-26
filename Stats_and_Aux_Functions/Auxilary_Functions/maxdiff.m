function maxval = maxdiff( X, Y )
% MAXDIFF( X, Y, p ) computes the maximal distance between two entries of
% the arrays X and Y.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   X   an array.
%   Y   an array of the same size as Y.
%
%--------------------------------------------------------------------------
% OUTPUT
%   maxval  the maximums norm of X or distance between X and Y.  
%--------------------------------------------------------------------------
% DEVELOPER TODOs:
% -------------------------------------------------------------------------
% EXAMPLES
% % % Define two arrays to compare
% X = [1, 2, 3, 4];
% Y = [5, 6, 7, 8];
% 
% % Find the maximum difference between elements of X and Y
% maxval = maxdiff(X, Y);
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow  
%--------------------------------------------------------------------------


%% ------------------------------------------------------------------------
%  check mandatory input
%--------------------------------------------------------------------------
if numel(X) ~= numel(Y)
   % default option of Y
   error( 'X and Y need to have the same number of elements.' )
end

%% ------------------------------------------------------------------------
%  main function
%--------------------------------------------------------------------------
maxval =  max( abs( X(:) - Y(:) ) );

return