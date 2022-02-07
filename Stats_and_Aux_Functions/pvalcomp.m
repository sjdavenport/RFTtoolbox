function [ pval_store ] = pvalcomp( dist1, dist2, doplot )
% NEWFUN serves as a function template.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% pvalcomp(normrnd(0,1,1, 10000), normrnd(0,1,1, 10000))
% pvalcomp(normrnd(0,1,1, 1000), chi2rnd(3,1, 10000))
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'doplot', 'var' )
   % Default value
   doplot = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
nvals = length(dist1);
pval_store = zeros(1, nvals);
% peak_height_sim_distbn = mean(local_maxima) - mean(peak_height_sim_distbn) + peak_height_sim_distbn;
for I = 1:nvals
    modul(I,1000)
    pval_store(I) = sum(dist2 > dist1(I));
end
pval_store = pval_store/length(dist2);

if doplot
    histogram(pval_store)
end

end

