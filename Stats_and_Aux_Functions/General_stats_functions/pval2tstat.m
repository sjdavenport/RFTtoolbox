function tstat = pval2tstat( pval, df )
% tstat2pval( tstat, do2sample ) calculates the pvalues from the
% t-statistic. 
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  tstat: an array of size dim
%  df: degrees of freedom
% Optional
%  do2sample: 0/1, 1: compute two sample p-values. 0: compute onesample
%                 p-values. Default is 1 ie to do two sample.
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
if ~exist( 'do2sample', 'var' )
   % Default value
   do2sample = 1;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
tstat = -tinv(pval, df);

end

