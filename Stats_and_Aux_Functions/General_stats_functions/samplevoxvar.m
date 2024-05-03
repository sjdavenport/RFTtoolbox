function out_var = samplevoxvar( array )
% SAMPLEVOXVAR( array ) finds the sample covariance of an array.
%--------------------------------------------------------------------------
% ARGUMENTS
% array     an nvox by D by nsubj array
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% mu = [0,0];
% Sigma = eye(2);
% data = mvnrnd(mu, Sigma, 100)';
% out = samplevoxcov(data)
% 
% nvox = 5;
% nsubj = 100;
% data = zeros(nvox,2,100);
% for I = 1:nvox
%     data(I,:,:) = mvnrnd(mu, Sigma, 100)'
% end
% out = samplevoxcov(data)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport

sA = size(array);

nsubj = sA(end);

mean_array = mean(array, length(sA));

out_var = sum((array - mean_array).*array, length(sA))/(nsubj - 1);

% nvox = sA(1);
% D = sA(2);
% out_cov = zeros(nvox, D^2);
% 
% for diagonal = 0:(D-1)
%     out_cov(:,diagonal
% end

end

