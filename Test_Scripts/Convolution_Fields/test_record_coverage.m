%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the record_coverage function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D examples
%% Simple 1D example
FWHM = 3; sample_size = 50; nvox = 5; mask = true(nvox, 1); resadd = 25;
spfn = @(nsubj) normrnd(0,1,[nvox, nsubj]); niters = 1000;
record_coverage( spfn, sample_size, FWHM, resadd, mask, niters)

%% %% 2D examples
%% small 2D example
FWHM = 3; sample_size = 50; Dim = [5,5]; mask = true(Dim); resadd = 25;
spfn = @(nsubj) normrnd(0,1,[Dim, nsubj]); niters = 1000;
record_coverage( spfn, sample_size, FWHM, resadd, mask, niters)

%% small 2D example - using HPE
FWHM = 3; sample_size = 50; Dim = [5,5]; mask = true(Dim); resadd = 11;
spfn = @(nsubj) normrnd(0,1,[Dim, nsubj]); niters = 1000;
record_coverage( spfn, sample_size, FWHM, resadd, mask, niters, 'hpe')

%% large 2D example
FWHM = 3; sample_size = 50; Dim = [50,50]; mask = true(Dim); resadd = 1;
spfn = @(nsubj) normrnd(0,1,[Dim, nsubj]); niters = 1000;
record_coverage( spfn, sample_size, FWHM, resadd, mask, niters)

%% 2D example with mask
 
%% %% 3D examples
%% Small 3D example
FWHM = 3; sample_size = 50; Dim = [5,5,5]; mask = true(Dim); resadd = 1;
spfn = @(nsubj) normrnd(0,1,[Dim, nsubj]); niters = 1000;
record_coverage( spfn, sample_size, FWHM, resadd, mask, niters)