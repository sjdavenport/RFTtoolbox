% test the 
addpath(genpath("/home/drtea/matlabToolboxes/RFTtoolbox/"))

clear all
close all

Dim = [30 30 30];
nSubj = 200;
FWHM = 4;
FWHMcor = 10;

% different choices for voxelmaps
rng(1)
voxelmap = 1:85184; % identity
voxelmap = randsample( 85184, 85184 ); % all random
voxelmap = (0:44:(44^3-1) ) + randsample( 44, 44 );
voxelmap = voxelmap(:);
voxelmap = (0:44:(44^3/2-1) ) + randsample( 44, 44 );
voxelmap = voxelmap(:);
voxelmap = [voxelmap; 44^3/2 + randsample( 44^3/2, 44^3/2 )];


rng(1)
[ data, RawNoise, TrnInd ] = noisegen_nonstat( Dim, nSubj, FWHM, FWHMcor, voxelmap );
figure(1)
imagesc(squeeze(data(:,14,:,14))), colorbar
figure(2)
imagesc(squeeze(RawNoise(:,14,:,14))), colorbar
figure(3)
imagesc(squeeze(data(:,14,:,15))), colorbar
figure(4)
imagesc(squeeze(RawNoise(:,14,:,15))), colorbar

% mean variance of nonstat field
datastd = std( data, 0, length(Dim)+1 );
figure(5), imagesc( squeeze(datastd(:,10,:)) ), colorbar
datamean = mean( data, length(Dim)+1 );
figure(6), imagesc( squeeze(datamean(:,10,:)) ), colorbar

% mean variance of raw data
rawstd = std( RawNoise, 0, length(Dim)+1 );
figure(7), imagesc( squeeze(rawstd(:,10,:)) ), colorbar
rawmean = mean( RawNoise, length(Dim)+1 );
figure(8), imagesc( squeeze(rawmean(:,10,:)) ), colorbar

% pointwise normalized fields
datan = (data - mean(data, 4) ) ./ std(data, 0, 4);

figure(9)
imagesc( squeeze(datan(:,14,:,14)) ), colorbar
figure(10)
imagesc(squeeze(datan(:,14,:,15))), colorbar
figure(11)
imagesc(squeeze(datan(:,14,:,1))), colorbar
figure(12)
imagesc(squeeze(datan(:,14,:,12))), colorbar