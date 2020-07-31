%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests to see what the difference is in masking or not
%%%    msaking the data on the LKC estimation
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D example


%% %% 2D Examples
%% Simple 2D example
%% small 2D sphere example
FWHM = 3; sample_size = 50; Dim = [5,5]; resadd = 3;
mask = bndry_voxels(true(Dim), 'full'); nsubj = 20;
lat_data = wnfield(mask, nsubj); 
params = ConvFieldParams( [FWHM,FWHM], resadd );

% Masking
cfield  = convfield( lat_data, params, 0 );
dcfield = convfield( lat_data, params, 1 );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );
L_conv

% Inconsistent masking
% This is wrong but is what I did once so wanted to see what effect is has
params.lat_masked = false;
cfield = convfield( lat_data, params, 0 );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );
L_conv
L0_conv %This is 0! so definitely wrong!

% No masking
params.lat_masked = false;
cfield = convfield( lat_data, params, 0 );
dcfield = convfield( lat_data, params, 1 );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );
L_conv
L0_conv

%% %% 3D Examples
%% MNImask
%% small 2D sphere example
FWHM = 3; resadd = 1;
mask = logical(imgload('MNImask')); nsubj = 20;
lat_data = wnfield(mask, nsubj); 
params = ConvFieldParams( [FWHM,FWHM,FWHM], resadd );

% Masking
cfield  = convfield( lat_data, params, 0 );
dcfield = convfield( lat_data, params, 1 );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );
L_conv

% Inconsistent masking
% This is wrong but is what I did once so wanted to see what effect is has
params.lat_masked = false;
cfield = convfield( lat_data, params, 0 );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );
L_conv
L0_conv %This is 0! so definitely wrong!

% No masking
params.lat_masked = false;
cfield = convfield( lat_data, params, 0 );
dcfield = convfield( lat_data, params, 1 );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );
L_conv
L0_conv