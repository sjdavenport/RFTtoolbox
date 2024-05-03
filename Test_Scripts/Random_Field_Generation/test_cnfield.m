%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the cnfield function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
% Simple 1D example
rng(99)
dim = [30 1];
voxmap = 1:dim;
% tmp2 = mod( voxmap, 6);
%voxmap = [ voxmap( tmp2 < 3 ), voxmap( tmp2 >= 3 ) ];
voxmap = randsample(dim(1), dim(1));

nsubj = 10;
FWHM  = 4; 
lat_data = cnfield( dim, 10, voxmap, 0, nsubj );

resadd = 3;
params = ConvFieldParams(FWHM, resadd);
smooth_data = convfield(lat_data, params);

subplot(1, 2, 1)
plot( lat_data(:,5) )

covest = rfcovest( smooth_data );

subplot(1, 2, 2)
imagesc(covest)

%% %% 2D Examples
dim    = [ 50 50 ];
nsubj  = 10
voxmap = 1:prod( dim )
% tmp2   = mod( voxmap, 12);
% voxmap = [ voxmap( tmp2 < 3 ), voxmap( tmp2 >= 3 & tmp2 < 6 ),...
%            voxmap( tmp2 >= 6 & tmp2 < 9 ),...
%            voxmap( tmp2 >= 9) ]
%    
cut   = 6;
shift = 4;
tmp2   = mod( voxmap, cut);
voxmap2 = [ voxmap( tmp2 <= cut-shift ), voxmap( tmp2 > cut-shift ) ]

lat_data = cnfield( dim, 10, voxmap2, 0, nsubj )
figure(1);
subplot(1,3,1)
imagesc( lat_data( :, :, 1 ) );
subplot(1,3,2)
imagesc( lat_data( :, :, 3 ) );
subplot(1,3,3)
imagesc( lat_data( :, :, 5 ) );


%% Better 2D example
%% %% 2D Examples
dim    = [ 50 50 ];
nsubj  = 10
voxmap = randsample(prod(dim),prod(dim),0)

lat_data = cnfield( dim, 10, voxmap, 3, nsubj )
figure(1);
subplot(1,3,1)
imagesc( lat_data( :, :, 1 ) );
subplot(1,3,2)
imagesc( lat_data( :, :, 3 ) );
subplot(1,3,3)
imagesc( lat_data( :, :, 5 ) );

