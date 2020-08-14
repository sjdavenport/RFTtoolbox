%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the EulerCharCrit function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D example
FWHM = 6; resadd = 1; nsubj = 50; nvox = 10;
lat_data = wnfield(nvox, 1);
a = EulerCharCrit( lat_data.field, lat_data.D, lat_data.mask );
a{1}
plot(lat_data);

%% Matlab version
a = EulerCharCrit( lat_data.field, lat_data.D, lat_data.mask, 'notC' );
a{1}

%% %% 2D Examples
%% Simple 2D example
FWHM = 6; resadd = 1; nsubj = 50; Dim = [50,50];
lat_data = wnfield(Dim, nsubj);
tic
a = EulerCharCrit( lat_data.field, lat_data.D, lat_data.mask );
a{1}
toc

%% tstat example
MNImask = imgload('MNImask');
MNImask_2D = logical(MNImask(:,:,45)); nsubj = 20;
lat_data = wnfield(MNImask_2D, nsubj);
params = ConvFieldParams( [3,3], 1);
tcfield = convfield_t( lat_data, params );
ECC = EulerCharCrit( tcfield.field, tcfield.D, tcfield.mask )
plot(ECC{1}(2:end-1,1), ECC{1}(3:end,2))
%%
imagesc(tcfield.field > -3)

%% Matlab version
tic
b = EulerCharCrit( lat_data.field, lat_data.D, lat_data.mask, 'notC' );
isequal(a,b)
toc

%% %% 3D Examples
%% Simple 2D example
FWHM = 6; resadd = 1; nsubj = 1; Dim = [70,90,70];
lat_data = wnfield(Dim, nsubj);
tic
b = EulerCharCrit( lat_data.field, lat_data.D, lat_data.mask );
toc

%% Matlab version
tic
a = EulerCharCrit( lat_data.field, lat_data.D, lat_data.mask, 'notC' );
isequal(a,b)
toc
