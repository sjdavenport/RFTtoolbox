%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the vectcov function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D examples
%% Simple 1D example
nsubj = 1000;
vect1 = normrnd(0,1,1,nsubj); vect2 = normrnd(0,1,1,nsubj);
vectcov( vect1, vect2, 2) %Should be about 0 

mu =[0,0]'; Sigma = [1,0.5;0.5,1];
data = mvnrnd(mu,Sigma,nsubj)';
vectcov(data(1,:), data(2,:) ) %Should be about 0.5

%% 
Sigma1 = [1,0.5;0.5,1]; Sigma2 = [1,3/4;3/4,1];
voxel1 = mvnrnd(mu,Sigma1,nsubj)'; voxel2 = mvnrnd(mu,Sigma2,nsubj)';
mate1 = zeros(2,nsubj); mate2 = mate1;
mate1(1,:) = voxel1(1,:); mate2(1,:) = voxel1(2,:);
mate1(2,:) = voxel2(1,:); mate2(2,:) = voxel2(2,:);
vectcov(mate1, mate2, 2)

%% %% 2D Examples
%% Simple 2D example
Dim = [5,5]; nsubj = 1000; mu =[0,0]'; Sigma = [1,0.5;0.5,1];
data = mvnrnd(mu,Sigma,nsubj*prod(Dim))';
mate1 = reshape(data(1,:), [Dim, nsubj]); mate2 = reshape(data(2,:), [Dim, nsubj]);
vectcov(mate1, mate2, 3)

%% More interesting 2D example
Sigma1 = [1,0.5;0.5,1]; Sigma2 = [1,3/4;3/4,1];
voxel1 = mvnrnd(mu,Sigma1,nsubj)'; voxel2 = mvnrnd(mu,Sigma2,nsubj)';
voxel3 = mvnrnd(mu,Sigma2,nsubj)'; voxel4 = mvnrnd(mu,Sigma1,nsubj)';
mate1 = zeros([2,2,nsubj]); mate2 = mate1;
mate1(1,1,:) = voxel1(1,:); mate2(1,1,:) = voxel1(2,:);
mate1(2,1,:) = voxel2(1,:); mate2(2,1,:) = voxel2(2,:);
mate1(1,2,:) = voxel3(1,:); mate2(1,2,:) = voxel3(2,:);
mate1(2,2,:) = voxel4(1,:); mate2(2,2,:) = voxel4(2,:);
vectcov(mate1, mate2)