function lap_rvs = rlap( b, mate_size, mu )
% rlap( b, mate_size, mu ) generates data from a laplacian distribution
% with scale parameter b and mean mu
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  b        the scale parameter
%  mate_size    the size of the matrix of random laplacian variable to
%               return
% Optional
%  mu       the mean of the distribution. Default is zero.
%--------------------------------------------------------------------------
% OUTPUT
% lap_rvs   a matrix of size mate_size of indepedent laplacian random
%           variables with scale b and mean mu
%--------------------------------------------------------------------------
% EXAMPLES
% lrvs = rlap( 1, [1,1000]);
% histogram(lrvs)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'mu', 'var' )
   % default option of opt1
   mu = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
X = exprnd(1/b, mate_size);
Y = exprnd(1/b, mate_size);
lap_rvs = X-Y + mu;

end

