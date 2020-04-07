% test the 
addpath(genpath("/home/drtea/matlabToolboxes/RFTtoolbox/"))

clear all
close all

Dim = [30 30 30];
nSubj = 50;
FWHM = 4;
FWHMcor = 5;

% different choices for voxelmaps
rng(1)
voxelmap = 1:85184; % identity
voxelmap = randsample( 85184, 85184 ); % all random
voxelmap = (0:44:(44^3-1) ) + randsample( 44, 44 );
voxelmap = voxelmap(:);

rng(1)
[ data, RawNoise, TrnInd ] = noisegen_nonstat( Dim, nSubj, FWHM, FWHMcor, voxelmap );
figure(1)
imagesc(data(:,:,14,14)), colorbar
figure(2)
imagesc(RawNoise(:,:,14,14)), colorbar
figure(3)
imagesc(data(:,:,14,15)), colorbar
figure(4)
imagesc(RawNoise(:,:,14,15)), colorbar

% mean variance of nonstat field
datastd = std( data, 0, length(Dim)+1 );
figure(5), imagesc(datastd(:,:,10)), colorbar
datamean = mean( data, length(Dim)+1 );
figure(6), imagesc( data(:,:,10)), colorbar

% mean variance of raw data
rawstd = std( RawNoise, 0, length(Dim)+1 );
figure(7), imagesc(rawstd(:,:,10)), colorbar
rawmean = mean( RawNoise, length(Dim)+1 );
figure(8), imagesc(rawmean(:,:,10)), colorbar