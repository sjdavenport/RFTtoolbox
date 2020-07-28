function threshold = spm_thresh( L, L0 )
% spm_thresh( L, L0 )
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

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'opt1', 'var' )
    % default option of opt1
    opt1 = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
% Convert the LKCs to resels
resel_vec = LKC2resel(L, L0);

% Calculate the RFT threshold using SPM
threshold = spm_uc_RF_mod(0.05,[1,sample_size-1],'T',resel_vec,1);

end

