%This file finds local maxima using spm_max
global stdsize
estmean = zeros(stdsize);
nSubj = 50;
for I = 1:nSubj
    estmean = estmean + readimg(I);
end
% surf(estmean(:,:,50))
estmean = estmean(:)/nSubj;
nVox = 91*109*91;

TmpMask = ones(stdsize); %Count all voxels for now!

[x,y,z] = ind2sub(stdsize,find(TmpMask));
XYZ     = [x y z]';

tic
[N,Z,M] = spm_max(estmean,XYZ);
toc

[~, sidx] = sort(Z, 'descend');
Z(sidx(1:20))
topMlocs = M(:,sidx(1:20));
lm_indices = sub2ind(stdsize,topMlocs(1,:),topMlocs(2,:),topMlocs(3,:));
estmean(lm_indices) %Gives the same as Z(sidx(1:20))
%So sub2ind(stdsize,topMlocs(1,:),topMlocs(2,:),topMlocs(3,:)) are the
%indices of estmean(:) (the vector rather than the matrix!) such that 
