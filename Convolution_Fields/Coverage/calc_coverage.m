function [voxcoverage, clustercoverage] = calc_coverage( maxnmin, LKCs, do2sample, alpha)
% CALC_COVERAGE( maxnmin, LKC_vec, alpha ) calculates the coverage that
% results from doing alpha level coverage. 
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
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
if ~exist( 'alpha', 'var' )
   % default option of opt1
   alpha = 0.05;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
niters = length(LKCs.L0);
thresholds = zeros(1, niters);

if do2sample
    alpha = alpha/2;
end

df = maxnmin.nsubj - 1;

for I = 1:niters
    thresholds(I) = EECthreshold( alpha, LKCs.L(:,I), LKCs.L0(I), "T", df );
end

peak_types = {'fine','lat','conv'};

% Obtain the voxelwise coverage
if do2sample
    for J = 1:3
        voxcoverage.(peak_types{J}) = sum(maxnmin.([peak_types{J},'maxima']) > thresholds);
        voxcoverage.(peak_types{J}) = voxcoverage.(peak_types{J}) + sum(maxnmin.([peak_types{J},'minima']) < -thresholds);
    end
else
    for J = 1:3
        voxcoverage.(peak_types{J}) = sum(maxnmin.([peak_types{J},'maxima']) > thresholds);
    end
end

% Obtain the clusterwise coverage
if do2sample
    for J = 1:3
        voxcoverage.(peak_types{J}) = sum(maxnmin.([peak_types{J},'maxima']) > thresholds);
        voxcoverage.(peak_types{J}) = voxcoverage.(peak_types{J}) + sum(maxnmin.([peak_types{J},'minima']) < -thresholds);
    end
else
    for J = 1:3
        voxcoverage.(peak_types{J}) = sum(maxnmin.([peak_types{J},'maxima']) > thresholds);
    end
end

end

