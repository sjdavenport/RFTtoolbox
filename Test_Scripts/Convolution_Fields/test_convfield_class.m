%%% 1D
%% Smoothing with increased resolution
nvox = 10; xvals = 1:nvox; FWHM = 2;
lat_data = normrnd(0,1,1,nvox); resadd = 10;
lat_field = fconv(lat_data, FWHM);
plot(xvals, lat_field, 'o-')
hold on
[convolution_field, xvals_fine] = convfield( lat_data, FWHM, resadd, 1);
plot(xvals_fine{1},convolution_field)

%% Multiple subjects
nsubj = 3; nvox = 100;
lat_data = normrnd(0,1,nvox,nsubj);
convolution_field = convfield( lat_data, FWHM, 0, 1 );
plot(1:nvox,convolution_field)

%% 1D derivatives
nvox = 100; resadd = 10; h = (1/(resadd+1)); 
lat_data = normrnd(0,1,nvox,1);
[convolution_field, xvals_fine] = convfield( lat_data, FWHM, resadd, 1);
deriv1 = convfield( lat_data, FWHM, resadd, 1, 1 );
deriv2 = diff(convolution_field)/h;
plot(xvals_fine{1}, deriv1)
hold on 
plot(xvals_fine{1}(1:end-1), deriv2, '--')

% 1D derivative (multiple subjects) (Some minor differences here!)
nvox = 100; FWHM = 3; xvals = 1:nvox; lat_data = normrnd(0,1,1,nvox);
aderiv = @(x) applyconvfield( x, lat_data, @(y) GkerMVderiv(y, FWHM) );
deriv_struct = convfield( lat_data, FWHM, resadd, 1, 1 );
deriv_conv = convfield( lat_data, FWHM, resadd, 1, 1 );
deriv_struct(1), deriv_conv(1), aderiv(1)

%% 2D
Dim = [25,25];
lat_data = normrnd(0,1,Dim);
smooth_data = convfield( lat_data, FWHM, 0, 2);
fine_data = convfield( lat_data, FWHM, 3, 2); %Convolution eval

zlimits = [min(fine_data(:))-0.1, max(fine_data(:))+0.1];

subplot(1,2,1)
surf(smooth_data)
zlim(zlimits)
title('Lattice Evaluation')
subplot(1,2,2)
surf(fine_data)
zlim(zlimits)
title('Convolution Field')

%% Dealing with the prefab issues
% smooth_data_prefab = convfield_prefab( lat_data, FWHM, 1, 2)
% smooth_data_pre = convfield_dep( lat_data, FWHM, 1, 2)

%% Matching to applyconvfield
cfield = @(x) applyconvfield(x, lat_data, FWHM);
mask = true(Dim);
cfieldnotrunc = @(x) applyconvfield(x, lat_data, FWHM, mask, 0);

smooth_data(20,20)
cfield([20,20]')
cfieldnotrunc([20,20]')
smooth_data(1,10)
cfield([1,10]')
cfieldnotrunc([1,10]')

%% 2D derivatives
Dim = [25,25];
lat_data = normrnd(0,1,Dim); resadd = 1;
derivfield = convfield( lat_data, FWHM, resadd, 2, 1);
surf(reshape(derivfield(:,:,1), spacep(Dim,resadd)))
title('2D 1st partial derivative of the convolution field')

%% Showing that the derivatives work
Dim = [5,5]; lat_data = normrnd(0,1,Dim);
point = [3,3]'; resadd = 100; h = 1/(1+resadd);

spaced_point = spacep(point,resadd);
derivfield = convfield( lat_data, FWHM, resadd, 2, 1);

% convolution derivatives
convfield_derivatives = squeeze(derivfield(spaced_point(1), spaced_point(2),:))

% Derivative using applyconvfield
aderiv = @(x) applyconvfield( x, lat_data, @(y) GkerMVderiv(y, FWHM)  );
acfield_derivatives = aderiv(point)

% Illustration on a fine lattice (not to be used in practice)
smoothfield100 = convfield( lat_data, FWHM, resadd, 2);
% smoothfield100 = convfield( lat_data, FWHM, resAdd, 2, 0);
partialderiv_finelat(1) = (smoothfield100(spaced_point(1)+1, spaced_point(2)) - smoothfield100(spaced_point(1),spaced_point(2)))/h;
partialderiv_finelat(2) = (smoothfield100(spaced_point(1), spaced_point(2) + 1) - smoothfield100(spaced_point(1),spaced_point(2)))/h;
fine_lattice_derivatives = partialderiv_finelat'
% note that it doesn't match perfectly because it's still a discrete
% approximation, but that's why we want to use derivfield in the first
% place!

% SPM (i.e. lattice) estimates of the derivative (quite off!)
smoothfield_spm = convfield( lat_data, FWHM, 0, 2);
spm_derivs(1) = (smoothfield_spm(point(1)+1, point(2)) - smoothfield_spm(point(1),point(2)));
spm_derivs(2) = (smoothfield_spm(point(1), point(2) + 1) - smoothfield_spm(point(1),point(2)));
spm_lattice_derivatives = spm_derivs'

%% 2D derivatives (multiple subjects)
Dim = [5,5]; nsubj = 20;
lat_data = normrnd(0,1,[Dim, nsubj]);
derivfield = convfield( lat_data, FWHM, 1, 2, 1)

%% 3D
%Compare to SPM
Dim = [10,10,10]; FWHM = 3;
lat_data = normrnd(0,1,Dim);
subplot(1,2,1)
spm_smooth_field = zeros(Dim); 
spm_smooth(lat_data, spm_smooth_field, FWHM)
surf(spm_smooth_field(:,:,5));
title('Lattice Eval')
subplot(1,2,2)
cfield = convfield( lat_data, FWHM, 0, 3); %Convolution eval
surf(cfield(:,:,5))
title('Convolution Field Eval (no smoothing)')

%% Fine evaluation
Dim = [10,10,10];
resadd = 10; D = length(Dim); FWHM = 3; 
slice = Dim(end)/2; spaced_slice = spacep(slice, resadd);
lat_data = normrnd(0,1,Dim);
cfield = convfield( lat_data, FWHM, resadd, D); %Convolution eval
twoDcfieldslice = cfield(:,:,spaced_slice);
zlimits = [min(twoDcfieldslice(:))-0.1, max(twoDcfieldslice(:))+0.1];

subplot(1,2,1)
spm_smooth_field = zeros(Dim); 
spm_smooth(lat_data, spm_smooth_field, FWHM)
surf(spm_smooth_field(:,:,slice));
title('Lattice Eval')
zlim(zlimits)

subplot(1,2,2)
surf(cfield(:,:,spaced_slice))
title('Convolution Field Eval (Convn)')
zlim(zlimits)

%% Compare to applyconvfield
lat_data = normrnd(0,1,Dim);
acfield = @(x) applyconvfield(x, lat_data, FWHM);
Dim = [10,10,10];
D = length(Dim); FWHM = 3; resadd = 0;
cfield = convfield( lat_data, FWHM, resadd, D); 
acfield([5,5,5]')
cfield(5,5,5)
acfield([1,1,10]')
cfield(1,1,10)

%% % 3D derivatives (1 subject)
Dim = [5,5,5]; D = length(Dim); FWHM = 3;
lat_data = normrnd(0,1,Dim); resadd = 18;
derivfield = convfield( lat_data, FWHM, 0, D, 1);
aderiv = @(x) applyconvfield( x, lat_data, @(y) GkerMVderiv(y, FWHM)  );
dfeval = squeeze(derivfield(3,3,3,:))
aceval = aderiv([3,3,3]')
spacing = 1/(1+resadd);
spaced_point = spacep( [3,3,3]', resadd);
cfield_fine = convfield( lat_data, FWHM, resadd, D);
pointeval = cfield_fine(spaced_point(1),spaced_point(2),spaced_point(3));
plusxeval = cfield_fine(spaced_point(1)+1,spaced_point(2),spaced_point(3));
plusyeval = cfield_fine(spaced_point(1),spaced_point(2)+1,spaced_point(3));
pluszeval = cfield_fine(spaced_point(1),spaced_point(2),spaced_point(3)+1);
derivx = (plusxeval - pointeval)/spacing;
derivy = (plusyeval - pointeval)/spacing;
derivz = (pluszeval - pointeval)/spacing;

fine_lat_deriv = [derivx,derivy,derivz]'

%% 3D derivatives (Multiple subjects)
Dim = [5,5,5]; D = length(Dim); FWHM = 3; nsubj = 2;
lat_data = normrnd(0,1,[Dim,nsubj]); resadd = 0;
derivfields = convfield( lat_data, FWHM, resadd, D, 1)
aderiv = @(x) applyconvfield( x, lat_data(:,:,:,2), @(y) GkerMVderiv(y, FWHM)  )
dfeval = squeeze(derivfields(3,3,3,2,:))
aceval = aderiv([3,3,3]')

%% Adjusting the field (1D)
nvox = 10; D = 1; FWHM = 2; lat_data = normrnd(0,1,1,nvox);
kernel = SepKernel( D, FWHM ); kernel.adjust = 0.1; resadd = 10;
[smoothfield, xvals_vecs] =  convfield( lat_data, FWHM, resadd, D);
[adjust_field, xvals_vecs_adjust] = convfield( lat_data, kernel, 0, D);

plot(xvals_vecs{1}, smoothfield)
hold on
plot(xvals_vecs_adjust{1},adjust_field)

acfield = @(tval) applyconvfield(tval, lat_data, FWHM);
adjust_field(3)
acfield(3.1)

%% Adjusting the field (3D)
Dim = [10,10,10]; D = length(Dim); FWHM = 1.5;
lat_data = normrnd(0,1,Dim); resadd = 9;
[smoothfield, xvals_vecs] = convfield( lat_data, FWHM, resadd, D );

kernel = SepKernel( D, FWHM ); kernel.adjust = [0.1,0,0];
[adjust_field, xvals_vecs_adjust] = convfield( lat_data, kernel, 0, D );

point = [1.1,1,1]'; spaced_point = spacep(point, resadd);
plot(xvals_vecs{1}, smoothfield(:,spaced_point(2), spaced_point(3)))
hold on
plot(xvals_vecs_adjust{1},adjust_field(:,point(2),point(3)))
acfield = @(tval) applyconvfield(tval, lat_data, FWHM);
adjust_field(1,1)
acfield(point)

%% %% Enlarging the field

%% 1D enlargement
nvox = 10; lat_data = normrnd(0,1,1,nvox); FWHM = 3; D = 1; resadd = 0;
cfield = convfield( lat_data', FWHM, resadd, D, 0, 1);
acfield = @(tval) applyconvfield(tval, lat_data, FWHM);

cfield(1)
acfield(0)

%% 1D enlargement
nvox = 10; lat_data = normrnd(0,1,1,nvox); FWHM = 3; D = 1; resadd = 1;
dx = 1/(1+resadd); enlarge = 1;
cfield = convfield( lat_data', FWHM, resadd, D, 0, enlarge);
acfield = @(tval) applyconvfield(tval, lat_data, FWHM);

cfield(1)
acfield(1-dx*enlarge)
%%
nvox = 10; lat_data = normrnd(0,1,1,nvox); FWHM = 3; D = 1; resadd = 0;
dx = 1/(1+resadd); enlarge = 5;
cfield = convfield( lat_data', FWHM, resadd, D, 0, enlarge);
acfield = @(tval) applyconvfield(tval, lat_data, FWHM);

cfield(1)
acfield(1-dx*enlarge)

%% 2D enlargement
Dim = [10,10]; lat_data = normrnd(0,1,Dim); FWHM = 3; D = 2; resadd = 1;
dx = 1/(1+resadd); enlarge = 1;
cfield = convfield( lat_data, FWHM, resadd, D, 0, enlarge)
acfield = @(tval) applyconvfield(tval, lat_data, FWHM);

cfield(1,1)
acfield([1-dx*enlarge,1-dx*enlarge]')

%% 3D enlargement
Dim = [10,10,10]; lat_data = normrnd(0,1,Dim); FWHM = 3; D = 3; resadd = 1;
dx = 1/(1+resadd); enlarge = 1;
cfield = convfield( lat_data, FWHM, resadd, D, 0, enlarge);
acfield = @(tval) applyconvfield(tval, lat_data, FWHM);

cfield(1,1,1)
acfield([1-dx*enlarge,1-dx*enlarge,1-dx*enlarge]')
