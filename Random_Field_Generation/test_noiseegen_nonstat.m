% test the 
addpath(genpath("/home/drtea/matlabToolboxes/RFTtoolbox/"))

clear all
close all

Dim = [30 30 30];
nSubj = 15;
FWHM = 2;
FWHMcor = 15;

% different choices for voxelmaps
rng(1)
Nwvox = 54872;
wd1 = round((Nwvox)^(1/3));

%voxelmap = 1:Nwvox; % identity
%voxelmap = randsample( Nwvox, Nwvox ); % all random
%voxelmap = (0:wd1:(Nwvox-1) ) + randsample( wd1, wd1 );
%voxelmap = voxelmap(:);
voxelmap = (0:wd1:(Nwvox/2-1) ) + randsample( wd1, wd1 );
voxelmap = voxelmap(:);
voxelmap = [voxelmap; Nwvox/2 + randsample( Nwvox/2, Nwvox/2 )];

%rng(1)
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

%% Test the multiplier base fields
data2 = multiplier_field( data, 1000 );
figure(13)
imagesc( squeeze(data2(:,14,:,14)) ), colorbar
figure(14)
imagesc(squeeze(data2(:,14,:,15))), colorbar
figure(15)
imagesc(squeeze(data2(:,14,:,1))), colorbar
figure(16)
imagesc(squeeze(data2(:,14,:,12))), colorbar

%% approximating isotropic fields by multiplier bootstrap
data3 = noisegen( Dim, 1000, FWHM );
data4 = multiplier_field( data3, 4 );

figure(17)
imagesc( squeeze( data3( :, 14, :, 1 ) ) ), colorbar
figure(18)
imagesc( squeeze( data4( :, 14, :, 2 ) ) ), colorbar
figure(19)
imagesc( squeeze( data4( :, 14, :,3 ) ) ), colorbar
figure(20)
imagesc( squeeze( data4( :, 14, :,4 ) ) ), colorbar
