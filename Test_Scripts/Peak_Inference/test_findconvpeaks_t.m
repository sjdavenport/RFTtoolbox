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

%% Small 2D application
Dim = [5,5]; nsubj = 50; lat_data = normrnd(0,1,[Dim,nsubj]); FWHM = 3;
resadd = 5; tfield = convfield_t( lat_data, FWHM, resadd );
surf(tfield)
findconvpeaks(lat_data, FWHM, 1, 'T')
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

%% %% Boundary examples
%% 2D
FWHM = 3; nsubj = 50; Dim = [5,5]; mask = true(Dim); 
lat_data = normrnd(0,1,[Dim, nsubj]); niters = 1000; resadd = 1; sresadd = 5;
[ tfield_lat ] = convfield_t( lat_data, FWHM, 0 );
% [ tfield_finelat, xvals ] = convfield_t( lat_data, FWHM, resadd );
[ tfield_superfinelat, xvals ] = convfield_t( lat_data, FWHM, 5 );

[ peak_est_locs, ~, peakvals ] = lmindices(tfield_lat, 3, mask);
[peaklox, peakvals] = findconvpeaks(lat_data, FWHM, peak_est_locs, 'T', mask)

max(tfield_superfinelat(:))
mask_hr = mask_highres(mask, sresadd);
[ finemaxlocs, ~, finemaxvals ] = lmindices(tfield_superfinelat, 3, mask_hr);
converted_locs = xvaleval(finemaxlocs, xvals);

%% Initialize on the fine lattice
[peaklox, peakvals] = findconvpeaks(lat_data, FWHM, converted_locs(:,1), 'T', mask)

%%
acf = @(tval) applyconvfield_t( tval, lat_data, FWHM, mask );
acf([5.5,0.5]')

lb = [0.5,0.5];
ub = [5.5,5.5];  

% There are no linear constraints, so set those arguments to |[]|. 
A = [];
b = [];
Aeq = [];
beq = [];  

% Choose an initial point satisfying all the constraints. 
x0 = peak_est_locs(:,1);  

% options = optimoptions(@fmincon,'Display','off'); % Ensures that no output 
                                                  % is displayed.
% options.Algorithm = 'sqp';
% mask_field( x, mask, xvals_vecs, 0 ) - 0.5;

mfield = @(x) mask_field(x, mask, 1:5, -1);

x = fmincon(@(y)-acf(y),x0,A,b,Aeq,beq,lb,ub,mfield)
%% 3D - Small example
FWHM = 3; nsubj = 50; Dim = [5,5,5]; mask = true(Dim); 
lat_data = normrnd(0,1,[Dim, nsubj]); niters = 1000; resadd = 1; sresadd = 9;
[ tfield_lat, xvals_lat ] = convfield_t( lat_data, FWHM, 1 );
% [ tfield_finelat, xvals ] = convfield_t( lat_data, FWHM, resadd );
[ tfield_superfinelat, xvals_super ] = convfield_t( lat_data, FWHM, sresadd );

[ peak_est_locs, ~, peakvals ] = lmindices(tfield_lat, 3, mask_highres(mask, 1));
peak_est_locs = xvaleval(peak_est_locs, xvals_lat);
[peaklox, peakvals] = findconvpeaks(lat_data, FWHM, peak_est_locs, 'T', mask)

max(tfield_superfinelat(:))
mask_hr = mask_highres(mask, sresadd);
[ finemaxlocs, ~, finemaxvals ] = lmindices(tfield_superfinelat, 3, mask_hr)
converted_locs = xvaleval(finemaxlocs, xvals_super);

%% 3D - Small example - with mask
FWHM = 3; nsubj = 50; Dim = [5,5,5]; mask = true(Dim); mask(2:4,2:4,2:4) = 0;
lat_data = normrnd(0,1,[Dim, nsubj]); niters = 1000; resadd = 1; sresadd = 1;
[ tfield_lat, xvals_lat ] = convfield_t( lat_data.*mask, FWHM, 1 );
% [ tfield_finelat, xvals ] = convfield_t( lat_data, FWHM, resadd );
[ tfield_superfinelat, xvals_super ] = convfield_t( lat_data.*mask, FWHM, sresadd );

[ peak_est_locs, ~, peakvals ] = lmindices(tfield_lat, 3, mask_highres(mask, 1));
peak_est_locs = xvaleval(peak_est_locs, xvals_lat);
[peaklox, peakvals] = findconvpeaks(lat_data, FWHM, peak_est_locs, 'T', mask)

max(tfield_superfinelat(:))
mask_hr = mask_highres(mask, sresadd);
[ finemaxlocs, ~, finemaxvals ] = lmindices(tfield_superfinelat, 3, mask_hr);
converted_locs = xvaleval(finemaxlocs, xvals_super)

%% 3D - Large example
FWHM = 3; nsubj = 50; Dim = [91,109,91]; mask = true(Dim); 
lat_data = normrnd(0,1,[Dim, nsubj]); niters = 1000; resadd = 1;
[ tfield_lat, xvals_lat ] = convfield_t( lat_data, FWHM, 1 );

[ peak_est_locs, ~, peakvals ] = lmindices(tfield_lat, 3, mask_highres(mask, 1));
peak_est_locs = xvaleval(peak_est_locs, xvals_lat);

tic
[peaklox, peakvals] = findconvpeaks(lat_data, FWHM, peak_est_locs, 'T', mask)
toc

tic
[peaklox, peakvals] = findconvpeaks_old(lat_data, FWHM, peak_est_locs, 'T', mask)
toc
%% Initialize on the fine lattice
[peaklox, peakvals] = findconvpeaks(lat_data, FWHM, converted_locs(:,1), 'T', mask)

%% Test to make sure it's finding the maxima
niters = 1000;
for I = 1:niters
    I
    FWHM = 3; nsubj = 50; Dim = [5,5,5]; mask = true(Dim);
    lat_data = normrnd(0,1,[Dim, nsubj]); niters = 1000; resadd = 1; sresadd = 9;
    [ tfield_lat, xvals_lat ] = convfield_t( lat_data, FWHM, 1 );
    % [ tfield_finelat, xvals ] = convfield_t( lat_data, FWHM, resadd );
    [ tfield_superfinelat, xvals_super ] = convfield_t( lat_data, FWHM, sresadd );
    
    [ peak_est_locs, ~, peakvals ] = lmindices(tfield_lat, 3, mask_highres(mask, 1));
    peak_est_locs = xvaleval(peak_est_locs, xvals_lat);
    [peaklox, peakvals] = findconvpeaks(lat_data, FWHM, peak_est_locs, 'T', mask);
    
    if max(tfield_superfinelat(:)) > max(peakvals) + 0.01
        disp('found one')
        
        disp('found one')
    end
end

%%
cfield = @(tval) applyconvfield_t( tval, lat_data, FWHM, mask)

cfield(peaklox(:,1))
cfield(peaklox(:,1)+[0.01,0,0]')
cfield(peaklox(:,1)+[0,0.01,0]')
cfield(peaklox(:,1)+[0,0,0.01]')

%%
mask_hr = mask_highres(mask, 5);
[ peak_est_locs, ~, peakvals ] = lmindices(tfield_superfinelat, 3, mask_hr);

%%
[peaklox, peakvals] = findconvpeaks_new(lat_data, FWHM, peak_est_locs, 'T', mask)
% xvaleval(peak_est_locs(:,1), 

%%
[peaklox, peakvals] = findconvpeaks_new(lat_data, FWHM, [4,1.5,2.5]', 'T', mask)

%% Simple 3D example with mask (implement to ensure it finds the peaks correctly)
Dim = [5,5,5]; nsubj = 25; FWHM = 3;
mask = true(Dim); mask(2:4,2:4,2:4) = 0;
lat_data = wnfield(mask, nsubj); resadd = 3;
[ tfield_finelat, cfields] = convfield_t_Field( lat_data, FWHM, resadd );

mask_hr = mask_highres(mask, resadd);
[peak_est_locs, ~,peakvals] = lmindices(tfield_finelat, 3, mask_hr);
peak_est_locs = xvaleval(peak_est_locs, cfields.xvals);
[peaklox, peakvals] = findconvpeaks(lat_data.field, FWHM, peak_est_locs, 'T', mask)

super_resadd = 11;
[cfield_fine, xvals_super] = convfield_t( lat_data.field.*mask, FWHM, super_resadd );
mask_superhr = mask_highres(mask, super_resadd);
[ finemaxlocs, ~, finemaxvals ] = lmindices(cfield_fine, 3, mask_superhr);
finemaxvals(1)
converted_locs = xvaleval(finemaxlocs, xvals_super);

%% Peak finding on the MNI mask
% Note you will need to have the MBS toolbox installed for this section to
% work
MNImask = logical(imgload('MNImask'));
nsubj = 25; FWHM = 3;
lat_data = wnfield(MNImask, nsubj); resadd = 1;
[ tfield_finelat, cfields] = convfield_t_Field( lat_data, FWHM, resadd );

mask_hr = mask_highres(MNImask, resadd);
[ finemaxlocs, ~, finemaxvals ] = lmindices(tfield_finelat, 1, mask_hr);
finemaxvals(1)
converted_locs = xvaleval(finemaxlocs, cfields.xvals);
[peaklox, peakvals] = findconvpeaks(lat_data.field, FWHM, converted_locs, 'T', MNImask)
