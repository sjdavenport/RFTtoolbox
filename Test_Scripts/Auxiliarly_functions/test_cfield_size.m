%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the size of convolution field
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1D 
nvox = 10; lat_data = randn([nvox,1]); D = 1; resadd = 0; FWHM = 3;
[cfield,xvals] = convfield( lat_data, FWHM, resadd, D, 0, ceil(resadd/2) );
voxel_locs = (ceil(resadd/2) + 1):(resadd+1):length(xvals{1})