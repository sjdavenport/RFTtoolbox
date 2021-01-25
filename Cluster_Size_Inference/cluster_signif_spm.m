function [ pval, EEC ] = cluster_signif_spm( cluster_extent, CDT, resel_vec )
% CLUSTER_SIGNIF_SPM( cluster_extent, CDT, resel_vec )
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% FWHM = 5; Dim = [250,250];
% resel_vec = spm_resels(FWHM, Dim, 'B')
% CDT = 2.5;
% pval = cluster_signif_spm( 50, CDT, resel_vec )
%
% FWHM = 5; Dim = [90,90,90];
% mask = ones(Dim); FWHM_vec = FWHM*ones(1, length(Dim));
% resel_vec =  spm_resels_vol(mask, FWHM_vec)';
% CDT = 2.5;
% [pval, EEC]  = cluster_signif_spm( 150, CDT, resel_vec )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
STAT = 'Z';
df = 1;

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'opt1', 'var' )
   % Default value
   opt1 = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
[pval, ~, EEC] = spm_P_RF(1,cluster_extent,CDT,df,STAT,resel_vec,1);

end

