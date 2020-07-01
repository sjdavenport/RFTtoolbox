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
%% 3D
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
%%
[peaklox, peakvals] = findconvpeaks_new(lat_data, FWHM, [4,1.5,2.5]', 'T', mask)

