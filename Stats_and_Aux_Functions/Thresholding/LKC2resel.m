function resel_vec = LKC2resel( LKCs, L_0 )
% LKC2RESEL( LKCs, L_0 ) takes in a vector of LKCs and returns the 
% corresponding vector of resels for use in SPM.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  LKCs     a vector of the LKCs such that L(i) is L_i i.e. the ith LKC for
%           i = 1:D. Note that L_0 is not included!
% Optional
%  L_0      the 0th LKC i.e. L_0. If this is not provided as an input it 
%           istaken to be 1
%--------------------------------------------------------------------------
% OUTPUT
% resel_vec     the vector of resels
%--------------------------------------------------------------------------
% DEVELOPER NOTES: Once we get rid of the SPM dependence this function will
% be deprecated
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Get rid of zero (high LKCs)
LKCs = LKCs(LKCs > 0);

% Ensure that the vector of LKCs is a row vector
if size(LKCs,1) > 1
    LKCs = LKCs';
end

% Compute the number of dimensions
D = length(LKCs);

%%  Main function
%--------------------------------------------------------------------------
% Calculate scaling factors to convert between resels and LKCs see e.g.
% Worsley 1992 (3D brain paper).
scaling_vec = repmat(sqrt(4*log(2)), 1, D).^(1:D);

% Initialise the resel vector
resel_vec = zeros(1,4);

% Set the non-zero LKCs
resel_vec(2:(D+1)) = LKCs./scaling_vec;

% Add in L_0 (defaulted to 1 if it has not been specified)
if ~exist('L_0', 'var')
    resel_vec(1) = L_0;
else
    resel_vec(1) = 1;
end

end

