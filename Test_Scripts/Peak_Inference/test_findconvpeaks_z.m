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

%% 1D with different xvals
Y = [1,1,2,2,1,1];
xvals = 11:(length(Y)+10);
Y = Field(Y', true(size(Y))'); Y.xvals{1} = xvals;
findconvpeaks(Y, 3, 1)

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
Y = [1,2,1,2,1]; FWHM = 2;
mask = logical([1,0,1,0,1]);
f = Field(Y', mask');
[maxloc, maxval] = findconvpeaks(f, FWHM, 1)

%% Raised the issues here with Fabian!
resadd = 10;
params = ConvFieldParams(4, resadd);
cfield  = convfield( f, params ); 
cfield = Mask(cfield);
plot(cfield.mask)

[maxloc, maxval] = findconvpeaks(f, FWHM, {1.5})
% Note that for a 1D convolution field the maximum can never lie on the
% boundary.

%% With the option for the mean!
Y = wfield(10,1)*(1/50);
x = 1:1:10; sigstore = peakgen1D( x, peakspec, peakparams, 1, smo );
plot(x, sigstore)
FWHM = 3; resadd = 1; params = ConvFieldParams(FWHM, resadd, 0);
meanonlat = convfield(sigstore, params);
truncation = 4*FWHM2sigma(FWHM);
meanfn = @(x) applyconvfield(x, sigstore, FWHM, true(1,length(sigstore)), truncation, meanonlat.xvals);

findconvpeaks(Y, 3, 1, 'Z', meanfn, meanonlat.field)
findlms(meanfn, 5, 4,6)

%% %% 2D Examples
%% Simple 2D application
Y = [1,1,1,1;1,2,2,1;1,2,2,1;1,1,1,1]; FWHM = 2; resadd = 50;
[maxloc, maxval] = findconvpeaks(Y, FWHM, 1)
params = ConvFieldParams([FWHM, FWHM], resadd);
fine_eval = convfield(Y, params);
max(fine_eval.field(:))
imagesc(fine_eval)

%% Finding peaks on the boundary
% An example where the peak will be outside the mask!
FWHM = 3; Y = [10,1,1;1,1,1;10,1,1]; 
mask = logical([1,1,1;0,1,1;1,1,1]);
f = Field(Y, mask);
[peaklocs, peakvals] = findconvpeaks(f, FWHM, [2,2]')

%%
% View masked field:
resadd = 21; params = ConvFieldParams([FWHM, FWHM], resadd);
cfield = Mask(convfield(Mask(f), params));
surf(cfield.field)
max(cfield.field(:))

applyconvfield([2,1.5]', f.field, FWHM, mask)

%% 2D multiple peaks
Y = [5,1,1,1;1,1,1,1;1,1,1,1;1,1,1,5]; resadd = 9;
FWHM = 2; params = ConvFieldParams([FWHM, FWHM], resadd);
cfield = convfield(Y, params);
surf(cfield.field)
[peaklocs, peakvals] = findconvpeaks(Y, FWHM, [1,1;4,4]');
% Note that the first column is the location of the first peak
% and the second column is the location of the second peak!

%% Finding peaks in the presence of a mask
Y = ones(4); FWHM = 2; resadd = 21;
params = ConvFieldParams([FWHM, FWHM], resadd);
mask =  logical([1,1,1,1;1,0,0,1;1,0,0,1;1,1,1,1]); 
f = Mask(Field(Y, mask));

field = @(tval) applyconvfield(tval, Y, FWHM, -1, 1:4, mask);

[maxloc, maxval] = findconvpeaks(f, FWHM, [1,1;1,4;4,1;4,4]')
cfield = Mask(convfield(f, params));
surf(cfield.field)
max(cfield.field(:))

%% %% 3D Examples
%% Simple 3D application (mask defaulted to ones everywhere!)
Y = ones([4,4,4]); Y(2:3,2:3,2:3) = 2; FWHM = 3;
[maxloc, maxval] = findconvpeaks(Y, FWHM, 1)

%% Large 3D application (takes a few seconds)
Dim = [91,109,91]; lat_data = randn(Dim); FWHM = 3;
resadd = 2; enlarge = ceil(resadd/2);
params = ConvFieldParams(repmat(FWHM, 1, 3), resadd, enlarge);

% Finds the maximum using findconvpeaks
[maxloc, maxval] = findconvpeaks(lat_data, FWHM) %Initializes on the integer lattice

% Evaluating on a lattice to illustrate that it does a good job!
[fine_field, xvals_vecs] = convfield(lat_data, params);
Dimhr = ( Dim - 1 ) * resadd + Dim + 2*enlarge';
[maxfineval, maxfineloc ] = max(fine_field.field(:));
converted_maxfineloc = convind(maxfineloc, Dimhr);
xvaleval(converted_maxfineloc', fine_field.xvals) 
% As there is a large image we have to take a low value of resadd here, 
% but as can be seen the true location is close to the peak on the resadd =
% 2 lattice.

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

