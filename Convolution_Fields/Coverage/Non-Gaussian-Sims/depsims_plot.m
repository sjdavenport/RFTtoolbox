global RFTboxloc
store_loc = [RFTboxloc, 'Convolution_Fields/Coverage/Non-Gaussian-Sims/'];

deep_coverage = load([store_loc, 'dep_coverage.mat']).coverage;
indeep_coverage = load([store_loc, 'indep_coverage.mat']).coverage;

nsubj_vec = 10:10:50;

nsims = min(size(deep_coverage,1), size(indeep_coverage,1)) - 1;
deep_coverage_store = zeros(nsims, length(nsubj_vec));
indeep_coverage_store = zeros(nsims, length(nsubj_vec));

for I = 1:nsims
    for J = 1:length(nsubj_vec)
        deep_coverage_store(I,J) = deep_coverage(I,J).coverage.conv;
        indeep_coverage_store(I,J) = indeep_coverage(I,J).coverage.conv;
    end
end

deep_std = std(deep_coverage_store);
indeep_std = std(indeep_coverage_store);
plot(nsubj_vec, deep_std);
hold on
plot(nsubj_vec, indeep_std);
legend('Dependent', 'Independent', 'Location', 'NW')
ylim([0, max([deep_std(:);indeep_std(:)])])
    
%%
nsubj = 40;
index = find(nsubj_vec == nsubj);
subplot(2,1,1)
h = histogram(deep_coverage_store(:,index), 'NumBins', 25);
binlims = h.BinLimits;
title('Dependent Coverage Range')
subplot(2,1,2)
histogram(indeep_coverage_store(:,index));
title('Independent Coverage Range')
xlim(binlims)
