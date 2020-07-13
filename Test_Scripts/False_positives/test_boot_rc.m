%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the boot_rc function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note that the bootstrap coverage values (for small nsubj) can be 
% considerably more variable than if they were calculated using
% independent data

%% %% 1D Examples
%% Simple 1D example
nvox = 50; nsubj = 100; FWHM = 3; data = randn(nvox, nsubj);
boot_rc( data, 10, FWHM, true(nvox,1) )

%% %% 2D Examples
%% Simple 2D example
Dim = [5,5]; nsubj = 100; FWHM = 3; sample_size = 10;
mask = true(Dim); data = randn([Dim, nsubj]);
boot_rc( data, sample_size, FWHM, mask )

%% %% 3D Examples
%% Simple 3D example
Dim = [5,5,5]; nsubj = 100; FWHM = 3; sample_size = 10;
mask = true(Dim); data = randn([Dim, nsubj]);
boot_rc( data, sample_size, FWHM, mask )

%% nifti 3D example
directory = ''; %(directory containing nifti files)
mask = imgload('MNImask'); %load MNI mask (requires the BrainSTAT toolbox)
sample_size = 10; resadd = 1; niters = 20;
boot_rc( directory, sample_size, FWHM, mask, resadd, niters )
