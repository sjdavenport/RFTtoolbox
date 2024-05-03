function threshold = spm_thresh( L, L0, field_type, df, alpha )
% spm_thresh( L, L0 ) calculates the RFT voxelwise threshold using SPM.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
%
%--------------------------------------------------------------------------
% EXAMPLES
%
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if ~exist('alpha', 'var')
    alpha = 0.05;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
% Convert the LKCs to resels
resel_vec = LKC2resel(L, L0);

% Calculate the RFT threshold using SPM
if strcmp(field_type, 'T')
    threshold = spm_uc_RF_mod(alpha,[1,df],'T',resel_vec,1);
else
    error('not implemented yet')
end

end

