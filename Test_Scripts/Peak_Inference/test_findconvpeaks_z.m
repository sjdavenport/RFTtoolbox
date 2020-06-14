%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the findconvpeaks function for mean 'Z' fields
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all 
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% Simple 1D application
Y = [1,2,1];
findconvpeaks(Y, 3, 1)

%% 1D with different xvals_vecs
Y = [1,1,2,2,1,1];
xvals_vecs = 11:(length(Y)+10);
findconvpeaks(Y, 3, 1, 'Z', ones(1,length(Y)), NaN, -1, xvals_vecs)

%% 1D multiple peaks
Y = [1,2,1,1,1,1,2,1];
findconvpeaks(Y, 3, 2) % Top 2 peaks
findconvpeaks(Y, 3, {1,6}) % Top 2 peaks starting initializing at 1, 6, 
%                                   needs the NaN in 1D to differentiate from 
%                                   the top number of peaks 
findconvpeaks(Y, 3, {4.5-0.00001}) %returns 4.5 which is not the max so watch
%                                   out it can sometimes get stuck, though
%                                   this is a very contrived example which
%                                   won't occur in reality

%% 1D with a mask
Y = [1,2,1,2,1]; FWHM = 2; D = 1;
mask = logical([1,0,1,0,1]);
[maxloc, maxval] = findconvpeaks(Y, FWHM, 1, 'Z', mask)
resadd = 5;
mask_hr = zero2nan(mask_highres(mask, resadd, ceil(resadd/2)));
cfield = mask_hr.*convfield( Y.*mask, FWHM, resadd, D, 0, ceil(resadd/2));
plot(cfield)

% Note that for a 1D convolution field the maximum can never lie on the
% boundary.

%% %% 2D Examples
%% Simple 2D application
Y = [1,1,1,1;1,2,2,1;1,2,2,1;1,1,1,1]; FWHM = 2; resadd = 10; D = 2;
[maxloc, maxval] = findconvpeaks(Y, FWHM, 1)
fine_eval = convfield(Y, FWHM, resadd, D);
max(fine_eval(:))
surf(fine_eval)

%% Finding peaks on the boundary
FWHM = 3; D = 2;
Y = [10,1,1;1,1,1;10,1,1]; %I.e. so the peak will be outside the mask!
mask = logical([1,1,1;0,1,1;1,1,1]);
[peaklocs, peakvals] = findconvpeaks(Y, FWHM, [2,2]', 'Z', mask)

% View masked field:
resadd = 9;
mask_hr = zero2nan(mask_highres(mask, resadd, ceil(resadd/2)));
cfield = mask_hr.*convfield( Y.*mask, FWHM, resadd, D, 0, ceil(resadd/2));
surf(cfield)
max(cfield(:))

%% 2D multiple peaks
Y = [5,1,1,1;1,1,1,1;1,1,1,1;1,1,1,5]; resadd = 9;
cfield = convfield(Y, 2, resadd, 2);
surf(cfield)
[peaklocs, peakvals] = findconvpeaks(Y, 2, [1,1;4,4]')

%% Finding peaks in the presence of a mask
Y = ones(4); FWHM = 2; D = 2;
mask =  logical([1,1,1,1;1,0,0,1;1,0,0,1;1,1,1,1]); resadd = 9;
field = @(tval) applyconvfield(tval, Y, FWHM, -1, 1:4, mask);
mask_hr = zero2nan(mask_highres(mask, resadd, ceil(resadd/2)));

[maxloc, maxval] = findconvpeaks(Y, FWHM, 1, 'Z', mask)
cfield = convfield(Y.*mask, FWHM, resadd, D, 0, ceil(resadd/2));
fine_eval = mask_hr.*convfield(Y.*mask, FWHM, resadd, D, 0, ceil(resadd/2));
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