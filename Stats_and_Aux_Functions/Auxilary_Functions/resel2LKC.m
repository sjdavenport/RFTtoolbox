function [ L, L0 ] = resel2LKC( resels )
% resel2LKC converts SPM resels to Lipshitz Killing curvatures
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  resels  a D + 1 length vector of resels
%--------------------------------------------------------------------------
% OUTPUT
%  L = [L1,..., L_D]
%  L0 (the zeroth LKC)
%--------------------------------------------------------------------------
% EXAMPLES
% FWHM = 3; Dim = [5,5];
% resels = spm_resels(FWHM,Dim, 'B')
% [L,L0] = resel2LKC(resels)
% resel_vec = LKC2resel(L,L0)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

% Ensure that resels is a row vector
if size(resels,2) == 1
    resels = resels';
end

%%  Main Function Loop
%--------------------------------------------------------------------------
D = find(resels,1,'last') - 1;
resels = resels(1:(D+1));

% Calculate scaling factors to convert between resels and LKCs see e.g.
% Worsley 1992 (3D brain paper).
scaling_vec = repmat(sqrt(4*log(2)), 1, D+1).^(0:D);

% Set the non-zero LKCs
scaled_resels = resels.*scaling_vec;
L = scaled_resels(2:(D+1));
L0 = scaled_resels(1);
end

