%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the mask_field function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note that the plots below are just included for illustration purposes, 
% mask_field is only designed for evaluating the high resolution mask at a 
% few points rather than everywhere (for that you should use mask_highres)

%% %% 1D example
nvox = 10; mask = true(1, nvox);
mask_field( [0, 0.5, 0.75, 3, nvox+0.5, nvox+1 ], mask)

%% %% 2D example
FWHM = 3; 
xvals_fine = 1:0.1:3; xvals_vecs = 1:3;
xvaluesatvoxels = xvals2voxels({xvals_fine,xvals_fine});

Y = [10,1,1;1,1,1;10,1,1]; %I.e. so the peak will be outside the mask!
mask = [1,1,1;0,1,1;1,1,1];

field = @(tval) applyconvfield(tval, Y, FWHM, mask );
mfield = @(x) mask_field(x, mask);
masked_field = @(x) mfield(x).*field(x);

fieldeval = reshape(masked_field(xvaluesatvoxels), [length(xvals_fine),length(xvals_fine)]);
surf(fieldeval)

%% Mask Viewing
% %% View masked field (old version)
% xvals_fine = 0.5:0.1:3.5;
% xvaluesatvoxels = xvals2voxels(xvals_fine,2)
% masked_field = @(x) mask_field(x, mask).*field(x);
% fieldeval = reshape(masked_field(xvaluesatvoxels), [length(xvals_fine),length(xvals_fine)]);
% surf(fieldeval)
% 
% Y = [1,1,1;10,1,1;1,1,1] %I.e so the peak will be outside the mask!
% mask = [0,1,1;0,1,1;0,1,1];
% findconvpeaks(Y, 3, [3,3]', mask) %Need to fix so that initialization at [2,2] is fine as well!
% 
% %% Display masked field using applyconvfield
% field = @(tval) applyconvfield(tval, Y, FWHM, -1, xvals, mask);
% masked_field = @(x) mask_field(x, mask).*field(x);
% xvals_fine = 0.5:0.1:(length(Y)+0.5);
% xvaluesatvoxels = xvals2voxels(xvals_fine,2);
% fieldeval = reshape(masked_field(xvaluesatvoxels), [length(xvals_fine),length(xvals_fine)]);
% surf(fieldeval)
% 
% %% View masked field (old version)
% FWHM = 3; D = 2;
% Y = [10,1,1;1,1,1;10,1,1]; %I.e. so the peak will be outside the mask!
% mask = logical([1,1,1;0,1,1;1,1,1]);
% xvals_fine = 0.5:0.1:3.5;
% xvaluesatvoxels = xvals2voxels(xvals_fine,2)
% field = @(tval) applyconvfield(tval, Y, FWHM, -1, 1:3, mask);
% 
% masked_field = @(x) mask_field(x, mask).*field(x);
% fieldeval = reshape(masked_field(xvaluesatvoxels), [length(xvals_fine),length(xvals_fine)]);
% surf(fieldeval)
% 
% Y = [1,1,1;10,1,1;1,1,1] %I.e so the peak will be outside the mask!
% mask = [0,1,1;0,1,1;0,1,1];
% findconvpeaks(Y, 3, [2,2]', mask) %Need to fix so that initialization at [2,2] is fine as well!
