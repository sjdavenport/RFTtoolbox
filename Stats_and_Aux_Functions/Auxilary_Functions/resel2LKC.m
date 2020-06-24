function LKCs = resel2LKC( resels )
% resel2LKC converts SPM resels to Lipshitz Killing curvatures
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  resels  a D + 1 length vector of resels
%--------------------------------------------------------------------------
% OUTPUT
% LKCs   a structure containing
%      LKCs.Lhat = L1,..., L_D
%      LKCs.L0 = L0 (the zeroth LKC)
%--------------------------------------------------------------------------
% EXAMPLES
% FWHM = 3; Dim = [5,5];
% resels = spm_resels(FWHM,Dim, 'B')
% LKCs = resel2LKC(resels)
% resel_vec = LKC2resel(LKCs)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Main Function Loop
%--------------------------------------------------------------------------
D = length(resels) - 1;

% Calculate scaling factors to convert between resels and LKCs see e.g.
% Worsley 1992 (3D brain paper).
scaling_vec = repmat(sqrt(4*log(2)), 1, D+1).^(0:D);

% Set the non-zero LKCs
scaled_resels = resels.*scaling_vec;
LKCs.hatL = scaled_resels(2:end);
LKCs.L0 = scaled_resels(1);
end

