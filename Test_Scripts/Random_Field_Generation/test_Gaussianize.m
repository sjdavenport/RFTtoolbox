%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the Gaussianize function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D example - tfield
lat_data = wtfield( [100,1], 1000, 3 ) 

mean(lat_data.field(:))
mean(lat_data.field(:))/std(lat_data.field(:))

gaussianized_data = Gaussianize( lat_data );
mean(gaussianized_data.field(:))
mean(gaussianized_data.field(:))/std(gaussianized_data.field(:))

subplot(2,1,1)
h = histogram(lat_data.field(1,:));
h.BinWidth = 0.25;
subplot(2,1,2)
h = histogram(gaussianized_data.field(1,:));
h.BinWidth = 0.25;

%% 1D example - original field is Gaussian
%% Simple 1D example
lat_data = wnfield( [100,1], 1000 ) 

mean(lat_data.field(:))
mean(lat_data.field(:))/std(lat_data.field(:))

gaussianized_data = Gaussianize( lat_data );
mean(gaussianized_data.field(:))
mean(gaussianized_data.field(:))/std(gaussianized_data.field(:))

subplot(2,1,1)
h = histogram(lat_data.field(1,:));
h.BinWidth = 0.25;
limits1 = h.BinLimits;
subplot(2,1,2)
h = histogram(gaussianized_data.field(1,:));
h.BinWidth = 0.25;
limits2 = h.BinLimits;
xlim([min(limits1(1),limits2(1)), max(limits1(2),limits2(2))])
    
%% %% 2D Examples
%% Simple 2D example
lat_data = wtfield( [100,100], 100, 3 );
tic
gaussianized_data = Gaussianize( lat_data );
toc
%%
vox = 1;
while 1
    vox = vox +1;
subplot(2,1,1)
histogram(lat_data.field(vox,1,:));
subplot(2,1,2)
h = histogram(gaussianized_data.field(vox,1,:));
h.BinWidth = 0.25;
pause
end
%% %% 3D Examples
%% Simple 3D example