%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the applyconvfield function
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
%% Simple 1D Example
Y = [1,2,3,4];
tval = 2; FWHM = 3;
applyconvfield(tval, Y, FWHM)

%% 1D Comparison with convfield
nvox = 50; lat_data = randn(1, nvox); FWHM = 3; resadd = 2;
acf = @(tval) applyconvfield(tval, lat_data, FWHM);
[cfield, xvals_vecs] = convfield( lat_data, FWHM, resadd, 1);

plot(xvals_vecs{1}, acf(xvals_vecs{1}))
hold on
plot(xvals_vecs{1}, cfield, '--');

%% 1D Comparison with spm\_conv (if installed)
nvox = 100; lat_field = normrnd(0,1,nvox,nsubj);
field_at_voxels = applyconvfield(1:nvox, mean(lat_field,2)', FWHM);
[~, ss] = applyconvfield(nvox/2, lat_field(:,1)', FWHM);
plot(field_at_voxels/sqrt(ss))
hold on
[smoothfield, ss_spm] = spm_conv_mod(mean(lat_field,2), FWHM);
plot(smoothfield'/sqrt(ss_spm));

%% 1D masking
FWHM = 4; nvox = 50; D = 1;
mask = true([1,nvox]); mask(20:40) = 0;
lat_data = randn(1,nvox); resadd = 2;
field = @(tval) applyconvfield( tval, lat_data, FWHM, mask );
mfield = @(x) mask_field(x, mask);
masked_field = @(x) mfield(x).*field(x);

[cfield, xvals_vecs] = convfield( lat_data, FWHM, resadd, D );

plot(xvals_vecs{1}, masked_field(xvals_vecs{1}))
hold on
plot(xvals_vecs{1}, cfield, '--');
legend('masked field', 'unmasked field')

%% %% 2D examples
%% Simple 2D example
Y = [1,2;3,4];
tval = [1.5,2,4; 3.4,2, 5]; FWHM = 3;
applyconvfield(tval, Y, FWHM)

%% Demonstrating the need for truncation (and for using applyconvfield)
% Need to truncate for speed, else  things are really slow!
FWHM = 3; lat_data = normrnd(0,1,1000,1000);
tic; applyconvfield([500,500]', lat_data, FWHM); toc
tic; applyconvfield([500,500]', lat_data, FWHM, 10); toc
tic; convfield(lat_data, FWHM, 0, 2); toc

%% %% 3D examples
%% Truncation or no truncation (that is the question)
noise = randn([91,109,91]); FWHM = 3; xvals_vecs = 1:91;
applyconvfield([50,40,60]', noise, FWHM, NaN, 0) % No truncation
sigma = FWHM2sigma(FWHM); truncation = round(4*sigma);
applyconvfield([50,40,60]', noise, FWHM, NaN, -1) %Default truncation
