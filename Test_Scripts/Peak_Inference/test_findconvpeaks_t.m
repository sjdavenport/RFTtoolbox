%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the findconvpeaks function for t fields
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all 
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D application
nvox = 10; nsubj = 20; lat_data = randn([nvox, nsubj]);
FWHM = 3; resadd = 3;
[ peakloc, peakval ] = findconvpeaks(lat_data, FWHM, 1, 'T')

[ tcfield, xvals_vecs ] = convfield_t( lat_data, FWHM, resadd );
plot(xvals_vecs{1}, tcfield)
xlim([xvals_vecs{1}(1), xvals_vecs{1}(end)])

tcf = @(tval) applyconvfield_t( tval, lat_data, FWHM );
tcf(peakloc -0.01)
tcf(peakloc +0.01)

%% 1D multiple peaks
nvox = 50; nsubj = 20; lat_data = randn([nvox, nsubj]); FWHM = 3; 
[ peakloc, peakval ] = findconvpeaks(lat_data, FWHM, 2, 'T')

resadd = 3; [ tcfield, xvals_vecs ] = convfield_t( lat_data, FWHM, resadd );
plot(xvals_vecs{1}, tcfield)
xlim([xvals_vecs{1}(1), xvals_vecs{1}(end)])

%% 1D with a mean
mu = [1,2,1,2,1]; FWHM = 2; nvox = length(mu); nsubj = 20; 
lat_data = mu' + randn([nvox, nsubj]);
[ maxloc, maxval ] = findconvpeaks(lat_data, FWHM, 1, 'T')
resadd = 20; [ tcfield, xvals_vecs ] = convfield_t( lat_data, FWHM, resadd );
plot(xvals_vecs{1}, tcfield)
xlim([xvals_vecs{1}(1), xvals_vecs{1}(end)])

%% 1D with a mask
FWHM = 5; nvox = length(mu); nsubj = 20;
lat_data = randn([nvox, nsubj]);
mu = [1,2,1,2,1]; lat_data = mu' + lat_data;
% lat_data = -mu' + lat_data + [1,2,1,100,1]';
mask = logical([1,0,1,0,1])';
[maxloc, maxval] = findconvpeaks(lat_data, FWHM, 1, 'T', mask)
resadd = 20;
mask_hr = zero2nan(mask_highres(mask, resadd, ceil(resadd/2)));
[ tcfield, xvals_vecs ] = convfield_t( lat_data.*mask, FWHM, resadd );
plot(xvals_vecs{1}, tcfield.*mask_hr)

%% %% 2D Examples
%% Simple 2D application
middle_value = 2; mu = zeros([4,4]) + 1; mu(2:3,2:3) = middle_value; FWHM = 2;
Dim = size(mu); nsubj = 20; lat_data = mu + randn([Dim, nsubj]);
[maxloc, maxval] = findconvpeaks(lat_data, FWHM, 1, 'T')
resadd = 10; fine_eval = convfield_t( lat_data, FWHM, resadd );
max(fine_eval(:))
surf(fine_eval)

%% Finding peaks on the boundary
FWHM = 3; corner_val = 10; mu = [corner_val,1,1;1,1,1;corner_val,1,1]; 
Dim = size(mu); nsubj = 20; lat_data = mu + randn([Dim, nsubj]);
mask = logical([1,1,1;0,1,1;1,1,1]);
[peaklocs, peakvals] = findconvpeaks(lat_data, FWHM, 1, 'T', mask)

% View masked field:
resadd = 5;
mask_hr = zero2nan(mask_highres(mask, resadd, ceil(resadd/2)));
tcfield = mask_hr.*convfield_t( lat_data.*mask, FWHM, resadd );
surf(tcfield)
max(tcfield(:))

%% 2D multiple peaks
mu = [5,1,1,1;1,1,1,1;1,1,1,1;1,1,1,5]; FWHM = 2;
Dim = size(mu); nsubj = 20; lat_data = mu + randn([Dim, nsubj]);
[peaklocs, peakvals] = findconvpeaks(lat_data, FWHM, 2, 'T')

resadd = 10; tcfield = convfield_t(lat_data, FWHM, resadd);
surf(tcfield)
max(tcfield(:))

%% Finding peaks in the presence of a central mask (need to improve the initialization ability!)
mu = ones(4); FWHM = 2;
Dim = size(mu); nsubj = 20; lat_data = mu + randn([Dim, nsubj]);
mask =  logical([1,1,1,1;1,0,0,1;1,0,0,1;1,1,1,1]); resadd = 9;
[ maxloc, maxval ] = findconvpeaks(lat_data, FWHM, 2, 'T', mask)

%%
resadd = 1;
tcfield = convfield_t( lat_data.*mask, FWHM, resadd );
mask_hr = zero2nan(mask_highres(mask, resadd, ceil(resadd/2)));
fine_eval = mask_hr.*tcfield;
surf(fine_eval)
max(fine_eval(:))

%% %% 3D Examples
%% Simple 3D application
Y = ones([4,4,4]); Y(2:3,2:3,2:3) = 2; FWHM = 3;
[maxloc, maxval] = findconvpeaks(Y, FWHM, 1)

%% Large 3D application
lat_data = randn([91,109,91]); FWHM = 3; D = 3; resadd = 2;
[maxloc, maxval] = findconvpeaks(lat_data, FWHM, 1)
fine_field = convfield(lat_data, FWHM, resadd, D, 0, ceil(resadd/2));
max(fine_field(:))

%% Finding peaks on the boundary
FWHM = 3; D = 3;
Y = ones([3,3,3]); Y(1,:,:) = 10; %I.e. so the peak will be outside the mask!
mask = true([3,3,3]); mask(1,2,2) = 0;
[peaklocs, peakvals] = findconvpeaks(Y, FWHM, [2,2,2]', 'Z', mask)

% View a slice through the masked field
resadd = 6;
mask_hr = zero2nan(mask_highres(mask, resadd, ceil(resadd/2)));
cfield = mask_hr.*convfield( Y.*mask, FWHM, resadd, D, 0, ceil(resadd/2));
surf(squeeze(cfield(1,:,:)))
max(cfield(:))