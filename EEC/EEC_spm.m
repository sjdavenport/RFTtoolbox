function EEC = EEC_spm( threshold, L, L0, field_type, nsubj )
% EEC_spm( LKCs, field_type ) calculates the expected Euler characteristic 
% of the excursion set above the given threshold
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  L    the postive LKCs
%  L0   the 0th LKC i.e. the Euler characteristic of the 
%  field_type  
%--------------------------------------------------------------------------
% OUTPUT
%  EEC   the expected Euler characteristic of the excursion set above the
%        given threshold
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

if ~exist('nsubj', 'var')
    nsubj = 1;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
% Uses SPM to calculate the expected Euler Characteristic (note the SPM
% functions are very confusing and we need to replace this but good to
% check against)
resel_vec = LKC2resel( L, L0 );

if strcmp(field_type, 'Z') 
    if ~isequal(nsubj, 1)
        error('For Gaussian fields need df = 1')
    end
elseif strcmp(field_type, 'T')
    df = [nsubj,1];
else
    error('fields other than Z and T are not coded')
end
[~,~,EEC] = spm_P_RF(1,0,threshold,df,field_type,resel_vec,1);

end

