global RFTboxloc
directory = [RFTboxloc, 'Convolution_Fields/Coverage/RS_testing/'];
even = load([directory, 'store_coverage_even.mat']);
odd = load([directory, 'store_coverage_odd.mat']);

nsubj_vec = 10:10:110;

conv_coverage_vec = zeros(1,length(nsubj_vec));
for I = 1:length(nsubj_vec)
   nsubj = nsubj_vec(I);
   if mod(nsubj,20) == 0 
       conv_coverage_vec(I) = even.store_coverage.(['nsubj_', num2str(nsubj)]).conv;
   else
       conv_coverage_vec(I) = odd.store_coverage.(['nsubj_', num2str(nsubj)]).conv;
   end    
end

plot(nsubj_vec, conv_coverage_vec)