%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the XXX function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D example (no padding)
nvox = 100; nsubj = 50; D = 1;
lat_data = wnfield(nvox, nsubj); FWHM = 2;
est_fwhm( lat_data, FWHM, 0 )
est_smooth(fconv(lat_data.field, FWHM, D))

%% Demonstrating the edge effect
nvox = 10; nsubj = 50; D = 1;
lat_data = wnfield(nvox, nsubj); FWHM = 6;
est_fwhm( lat_data, FWHM, 0 )
est_smooth(fconv(lat_data.field, FWHM, D))

%% 1D with padding
nvox = 10; nsubj = 50; mask = true(nvox,1); 
FWHM = 6; pad = ceil( 4 * FWHM2sigma( FWHM ) );
padded_mask = logical( pad_vals( mask, pad) );
lat_data = wnfield(size(padded_mask,1), nsubj); lat_data.mask = padded_mask;
est_fwhm( lat_data, FWHM, 0 )
est_smooth(fconv(lat_data.field, FWHM, D))

%% %% 2D Examples
%% Simple 2D example
Dim = [50,50]; nsubj = 50; D = length(Dim);
lat_data = wnfield(Dim, nsubj); FWHM = 2;
est_fwhm( lat_data, FWHM, 0 )
est_smooth(fconv(lat_data.field, FWHM, D))

%% Demonstrating the edge effect
Dim = [10,10]; nsubj = 50; D = length(Dim);
lat_data = wnfield(Dim, nsubj); FWHM = 20;
est_fwhm( lat_data, FWHM, 0 )
est_smooth(fconv(lat_data.field, FWHM, D))

%% 2D with padding (Once you get to FWHM = 20, it doesn't matter which method you use)
Dim = [10,10]; nsubj = 50; D = length(Dim); mask = true(Dim);
FWHM = 20; pad = ceil( 4 * FWHM2sigma( FWHM ) );
padded_mask = logical( pad_vals( mask, pad) );
lat_data = wnfield(padded_mask, nsubj); 
est_fwhm( lat_data, FWHM, 0 )
est_smooth(fconv(lat_data.field, FWHM, D), padded_mask)

%% %% 3D Examples
%% Simple 3D example (clear edge effect!)
Dim = [25,25,25]; nsubj = 50; D = length(Dim);
lat_data = wnfield(Dim, nsubj); FWHM = 3;
est_fwhm( lat_data, FWHM, 0 )
est_smooth(fconv(lat_data.field, FWHM, D))

%% 3D with padding (at this large FWHM, it doesn't matter what method you use)
Dim = [10,10,10]; nsubj = 50; D = length(Dim); mask = true(Dim);
FWHM = 10; pad = ceil( 4 * FWHM2sigma( FWHM ) );
padded_mask = logical( pad_vals( mask, pad) );
lat_data = wnfield(padded_mask, nsubj); 
est_fwhm( lat_data, FWHM, 0 )
est_smooth(fconv(lat_data.field, FWHM, D), padded_mask)
