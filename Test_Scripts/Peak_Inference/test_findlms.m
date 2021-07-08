%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the findlms function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all 
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D examples
%% Simple 1D example
FWHM = 3;
Y = [1,2,1];
findconvpeaks(Y, FWHM, 1)
cfield = @(tval) applyconvfield(tval, Y, 3)
findlms( cfield, 2.5, 1, 3 )

%% Multiple peaks - same height
Y = [1,2,1,1,1,1,1,2,1]; FWMM = 3;
findconvpeaks(Y, FWHM, 2)
cfield = @(tval) applyconvfield(tval, Y, FWHM)
% xvals_fine = 1:0.1:length(Y);
% plot(xvals_fine, convfield(Y, FWHM, 0.1, 1))
findlms( cfield, [2.5,6.5], {1,3}, {6,8})

%% %% 2D examples
%% Simple 2D example
FWHM = 3; Y = [1,1,1,1;1,2,2,1;1,2,2,1;1,1,1,1];
findconvpeaks(Y, FWHM, 1)

cfield = @(tval) applyconvfield(tval, Y, FWHM);
findlms( cfield, [2,2]', {[1.5,1.5]'}, {[3,3]'})

%% 2D multiple peaks
Y = [5,1,1,1;1,1,1,1;1,1,1,1;1,1,1,5]
findconvpeaks(Y, 2, [1,1;4,4]')
cfield = @(tval) applyconvfield(tval, Y, 2)
findlms( cfield, [1,1;4,4]', {[1,1]', [3,3]'},  {[2,2]', [4,4]'})

%%
% Works with functions that take NaN values so long as the initial
% estimate is well defined!
mask = [0,1,1];
mask2 = [0,1,0];
Y = 3:-1:1;
FWHM = 4;
mfield = @(x) zero2nan(mask_field( x, mask2 ));
cfield = @(x) applyconvfield(x, Y, FWHM, -1, 1:3, mask)
masked_field = @(x) -mfield(x).*cfield(x);
xvals_fine = 1:0.1:3;
plot(xvals_fine, masked_field(xvals_fine))
xlim([1,3])
findlms( masked_field, 2, 2 )

%% NaN example
mask = [0,1,1,0,1,1];
mask2 = [0,1,0,0,1,0];
Y = [3:-1:1, 3:-1:1];
nvox = length(Y);
FWHM = 4;
mfield = @(x) zero2nan(mask_field( x, mask2 ));
cfield = @(x) applyconvfield(x, Y, FWHM, mask); %Note we haven't multiplied Y by the mask here!
masked_field = @(x) -mfield(x).*cfield(x);
xvals_fine = 1:0.1:nvox;
plot(xvals_fine, masked_field(xvals_fine))
xlim([1,nvox])
findlms( masked_field, 2, 2 )
findlms( masked_field, 5, 10 ) % Note that it fails to search beyond the NaNs!

%% 2D example comparing box size and search time
Dim = [5,5]; nsubj = 50; lat_data = normrnd(0,1,[Dim,nsubj]); FWHM = 3;
tcfield = convfield_t(lat_data, FWHM, 0);
lmloc = lmindices(tcfield);
tcf = @(tval) applyconvfield_t( tval, lat_data, FWHM );
niters = 1000;
store_time = zeros(1,niters);
for I = 1:niters
    I
    tic
    findlms(tcf, lmloc, 1)
    store_time(I) = toc;
end
mean(store_time(I))
%0.25 for boxsize = 5
%0.29 for boxsize 1!