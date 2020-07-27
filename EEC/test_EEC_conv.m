% test the LKC convoultion estimate using Sam's idea of peaks between
% voxels from convolution fields
addpath(genpath("/home/drtea/matlabToolboxes/RFTtoolbox/"))

clear all
close all

Dim = [30 30 30];
D = length(Dim);
nSubj = 15;
FWHM = 2;
[ data, RawNoise ] = noisegen( Dim, nSubj, 5 );

figure(1)
imagesc(RawNoise(:,:,1,1)), colorbar
figure(2)
imagesc(data(:,:,1,1)), colorbar

%% %%%%%%%%%%%%%%%
EC = EulerCharCrit( data, D );
n = 1

% find the critical value locations from EC
critValsn = EC{n}(2:end-1,1);
critsLoc = [];

for i = 1:length(critValsn )
    ind = find(data(:,:,:,1) == critValsn(i));

    [I1,I2,I3] = ind2sub(Dim,ind);
    critsLoc = [ critsLoc, [I1;I2;I3]];    
end

peak_locs = findconvpeaks( RawNoise(:,:,:,n), FWHM, critsLoc(:,1:5) );