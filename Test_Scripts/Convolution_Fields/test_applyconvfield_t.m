%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the applyconvfield_t function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note that the plots below are just included for illustration purposes,
% applyconvfield_t is only designed for evaluating the field at a few
% points rather than everywhere

%% %% 1D examples
%% 1D t convolution field
FWHM = 4; nvox = 100; nsubj = 40;
lat_data = normrnd(0,1, nvox, nsubj); resadd = 2;
tcf = @(tval) applyconvfield_t( tval, lat_data, FWHM );
[cfield, xvals_vecs] = convfield_t( lat_data, FWHM, resadd );

plot(xvals_vecs{1}, tcf(xvals_vecs{1}))
hold on
plot(xvals_vecs{1}, cfield, '--');

%% 1D masked example
FWHM = 4; nvox = 50; nsubj = 20;
mask = true([1,nvox]); mask(20:40) = 0;
lat_data = normrnd(0,1, nvox, nsubj); resadd = 2;
tcf = @(tval) applyconvfield_t( tval, lat_data, FWHM, mask );

mfield = @(x) mask_field(x, mask);
masked_field = @(x) mfield(x).*tcf(x);

[cfield, xvals_vecs] = convfield_t( lat_data, FWHM, resadd );

plot(xvals_vecs{1}, masked_field(xvals_vecs{1}))
hold on
plot(xvals_vecs{1}, cfield, '--');
legend('masked field', 'unmasked field')

%% %% 2D examples (need to work on this)
% nsubj = 20;
% Dim = [10,10];
% xvals_vecs = {1:Dim(1), 1:Dim(2)};
% xvaluesatvoxels = xvals2voxels(xvals_vecs);
% lat_field = normrnd(0,1,[Dim, nsubj]);
% field_at_voxels = reshape(applyconvfield_t( xvaluesatvoxels, lat_field, FWHM, -1 ), Dim);
% smoothfield = smoothtstat( lat_field, FWHM );
% surf(field_at_voxels)
% pause
% surf(smoothfield)
% 
% nvox = 100; nsubj = 50; FWHM = 3;
% lat_data = normrnd(0,1,nvox,nsubj);
% findconvpeaks_t(lat_data, FWHM)
% xvals_fine = 1:0.1:nvox;
% tcf = @(tval) applyconvfield_t( tval, lat_data, FWHM );
% plot(xvals_fine, tcf(xvals_fine))
