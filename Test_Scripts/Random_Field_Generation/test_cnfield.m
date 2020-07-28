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
%% Simple 1D example
dim = [50 1];
voxmap = 1:dim
tmp2 = mod( voxmap, 6);
voxmap = [ voxmap( tmp2 < 2 ), voxmap( tmp2 >= 2 ) ]

lat_data = cnfield( dim, 10, voxmap, 0, 10 )
plot( lat_data(:,5) );

%% %% 2D Examples
dim    = [ 50 50 ];
nsubj  = 10
voxmap = 1:prod( dim )
tmp2   = mod( voxmap, 9);
voxmap = [ voxmap( tmp2 < 3 ), voxmap( tmp2 >= 3 & tmp2 < 6 ),...
           voxmap( tmp2 >= 6) ]

lat_data = cnfield( dim, 10, voxmap, 2, nsubj )
figure(1);
subplot(1,3,1)
imagesc( lat_data( :, :, 1 ) );
subplot(1,3,2)
imagesc( lat_data( :, :, 3 ) );
subplot(1,3,3)
imagesc( lat_data( :, :, 5 ) );

