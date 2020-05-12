function cov_array = vectcov( array1, array2, dimension, normalize )
% VECTCOV( array1, array2, dimension, normalize ) calculates the covariance
% of two arrays along a given dimension.
%--------------------------------------------------------------------------
% ARGUMENTS
% array1
% array2
% dimension
% normalize
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% % 1D 
% nsubj = 1000;
% vect1 = normrnd(0,1,1,nsubj); vect2 = normrnd(0,1,1,nsubj);
% vectcov( vect1, vect2, 2)
%
% mu =[0,0]'; Sigma = [1,0.5;0.5,1];
% data = mvnrnd(mu,Sigma,nsubj)';
% vectcov(data(1,:), data(2,:) )
%
% Sigma1 = [1,0.5;0.5,1]; Sigma2 = [1,3/4;3/4,1];
% voxel1 = mvnrnd(mu,Sigma1,nsubj)'; voxel2 = mvnrnd(mu,Sigma2,nsubj)';
% mate1 = zeros(2,nsubj); mate2 = mate1;
% mate1(1,:) = voxel1(1,:); mate2(1,:) = voxel1(2,:);
% mate1(2,:) = voxel2(1,:); mate2(2,:) = voxel2(2,:);
% vectcov(mate1, mate2, 2)
%
% % 2D
% Dim = [5,5]; nsubj = 1000; mu =[0,0]'; Sigma = [1,0.5;0.5,1];
% data = mvnrnd(mu,Sigma,nsubj*prod(Dim))';
% mate1 = reshape(data(1,:), [Dim, nsubj]); mate2 = reshape(data(2,:), [Dim, nsubj]);
% vectcov(mate1, mate2, 3)
% 
% Sigma1 = [1,0.5;0.5,1]; Sigma2 = [1,3/4;3/4,1];
% voxel1 = mvnrnd(mu,Sigma1,nsubj)'; voxel2 = mvnrnd(mu,Sigma2,nsubj)';
% voxel3 = mvnrnd(mu,Sigma2,nsubj)'; voxel4 = mvnrnd(mu,Sigma1,nsubj)';
% mate1 = zeros([2,2,nsubj]); mate2 = mate1;
% mate1(1,1,:) = voxel1(1,:); mate2(1,1,:) = voxel1(2,:);
% mate1(2,1,:) = voxel2(1,:); mate2(2,1,:) = voxel2(2,:);
% mate1(1,2,:) = voxel3(1,:); mate2(1,2,:) = voxel3(2,:);
% mate1(2,2,:) = voxel4(1,:); mate2(2,2,:) = voxel4(2,:);
% vectcov(mate1, mate2)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
s1 = size(array1); s2 = size(array2);
if nargin < 3
    dimension = length(s1);
end
if nargin < 4
    normalize = 1;
end

if ~isequal(s1, s2)
   error('The arrays must have the same dimension')
end
if dimension > length(s1)
    error('The dimension must be <= the number of dimensions of the arrays')
end

array1_mean = mean(array1, dimension);

cov_array = mean((array1 - array1_mean).*array2, dimension);

N = size(array1, dimension);

if normalize == 1
    cov_array = (cov_array)*(N/(N-1));
end

end
