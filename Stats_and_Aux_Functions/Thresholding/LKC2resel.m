function resel_vec = LKC2resel( L, L0 )
% LKC2RESEL( LKCs ) takes in a vector of LKCs and returns the 
% corresponding vector of resels for use in SPM.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  L = [L1,..., L_D]
%  L0 (the zeroth LKC)
%--------------------------------------------------------------------------
% OUTPUT
%  resel_vec     the vector of resels
%--------------------------------------------------------------------------
% DEVELOPER NOTES: Once we get rid of the SPM dependence this function will
% be deprecated
%--------------------------------------------------------------------------
% EXAMPLES
% FWHM = 3; Dim = [5,5];
% resels = spm_resels(FWHM,Dim, 'B')
% [L,L0] = resel2LKC(resels)
% resel_vec = LKC2resel(L,L0)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Get rid of zero (high LKCs)
if L(end) == 0
    resel_vec = LKC2resel( L(1:end-1), L0);
    return
end

% Ensure that the vector of LKCs is a row vector
if size(L,1) > 1
    L = L';
end

% Compute the number of dimensions
D = length(L);

%%  Main function
%--------------------------------------------------------------------------
% Calculate scaling factors to convert between resels and LKCs see e.g.
% Worsley 1992 (3D brain paper).
scaling_vec = repmat(sqrt(4*log(2)), 1, D).^(1:D);

% Initialise the resel vector
resel_vec = zeros(1,D+1);

% Set the non-zero LKCs
resel_vec(2:(D+1)) = L./scaling_vec;
resel_vec(1) = L0;
end

