function [voxcoverage, clustercoverage, thresholds] = calc_coverage( maxnmin, LKCs, do2tail, alpha, option)
% CALC_COVERAGE( maxnmin, LKC_vec, alpha ) calculates the coverage that
% results from doing alpha level coverage. 
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  maxnmin   a structure containing the maxima and minima information. This 
%            is obtained by running record_coverage
% Optional
%  do2tail   0/1 whether to do a one sample or a two sample test. Default
%              is 1, i.e. to do a 2 sample test
%  alpha     the alpha level at which to threshold. Default is 0.05.
%--------------------------------------------------------------------------
% OUTPUT
%  voxcoverage    a structure giving the voxelwise coverage on the lattice,
%                 fine lattice and of the convolution field
%  clustercoverage  a structure giving the cluster coverage on the lattice,
%                 fine lattice and of the convolution field
%--------------------------------------------------------------------------
% EXAMPLES
% [voxcoverage, clustercoverage] = calc_coverage( maxnmin, LKCs, do2sample)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'do2tail', 'var' )
   % default option of opt1
   do2tail = 1;
end

if ~exist( 'alpha', 'var' )
   % default option of opt1
   alpha = 0.05;
end

if ~exist('option', 'var')
    option = 'standard';
end

%%  Main Function Loop
%--------------------------------------------------------------------------
niters = min(length(maxnmin.alphathresholds),  size(LKCs.L,2));
thresholds = zeros(1, niters);

if do2tail
    alpha = alpha/2;
end

df = maxnmin.nsubj - 1;

for I = 1:niters
    thresholds(I) = EECthreshold( alpha, LKCs.L(:,I)', LKCs.L0(I), "T", df, option );
end

peak_types = {'finelat','lat','conv'};

% Obtain the voxelwise coverage
for J = 1:3
    setofmaxima = maxnmin.([peak_types{J},'maxima']);
    voxcoverage.(peak_types{J}) = sum(setofmaxima(1:niters) > thresholds);
end
if do2tail
    for J = 1:3
        setofminima = maxnmin.([peak_types{J},'minima']);
        voxcoverage.(peak_types{J}) = voxcoverage.(peak_types{J}) + sum(setofminima(1:niters) < -thresholds);
    end
end

% Scale by the number of iterations
for J = 1:3
    voxcoverage.(peak_types{J}) = voxcoverage.(peak_types{J})/niters;
end

% Obtain the clusterwise coverage
largerthanthresh = maxnmin.allmaxima(:,1:niters) > thresholds;
clustercoverage = sum(largerthanthresh(:));
if do2tail
    lowerthanthresh = maxnmin.allminima(:,1:niters) < -thresholds;
    clustercoverage = clustercoverage + sum(lowerthanthresh(:));
end

% Scale by the number of iterations
clustercoverage = clustercoverage/niters;

end

