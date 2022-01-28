%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the convfield_t function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1D convolution t field
nvox = 10; nsubj = 20; resadd = 20; FWHM = 2;
lat_data = normrnd(0,1,[nvox,nsubj]);
params = ConvFieldParams(FWHM, resadd);
tcfield = convfield_t(lat_data, params);
lattice_tfield = convfield_t(lat_data, FWHM);
plot(1:nvox, lattice_tfield.field, 'o-')
hold on
plot(tcfield.xvals{1}, tcfield.field)
title('1D convolution t field')
legend('Convolution field', 'Lattice Evaluation', 'Location', 'Best')
xlabel('voxels')
ylabel('t field')

%% 2D convolution field
Dim = [10,10]; nsubj = 20; resadd = 5; FWHM = 2;
params = ConvFieldParams([FWHM, FWHM], resadd);
lat_data = normrnd(0,1,[Dim,nsubj]);
tcfield = convfield_t(lat_data, params);
lattice_tfield = convfield_t(lat_data, FWHM);
zmax = max(tcfield.field(:)); zmin = min(tcfield.field(:));
subplot(1,2,1); imagesc(lattice_tfield); zlim([zmin,zmax])
title('2D lattice t field')
xlabel('x'); ylabel('y'); zlabel('t field')
subplot(1,2,2); imagesc(tcfield); zlim([zmin,zmax])
title('2D convolution t field')
xlabel('x'); ylabel('y'); zlabel('t field')

%% Matching to applyconvfield_t
D = 2; dim = [25,25]; mask = true( dim ); FWHM = 3; nsubj = 50;
params = ConvFieldParams( FWHM * ones( [ 1 D ] ), 0, 0, false );
params_fine = ConvFieldParams( FWHM * ones( [ 1 D ] ), 3, 0, false );
lat_data = wfield( mask, nsubj );
smooth_tfield = convfield_t( lat_data, params );
fine_tfield = convfield( lat_data, params_fine ); %Convolution eval

zlimits = [min(fine_data.field(:))-0.1, max(fine_data.field(:))+0.1];

subplot(1,2,1)
surf(smooth_data.field)
zlim(zlimits)
title('Lattice tfield Evaluation')
subplot(1,2,2)
surf(fine_data.field)
zlim(zlimits)
title('Convolution t Field')

%% Matching to applyconvfield
tcfield = @(x) applyconvfield_t(x, lat_data.field, FWHM);
mask = true(dim);
tcfieldnotrunc = @(x) applyconvfield_t(x, lat_data.field, FWHM, mask, 0);

smooth_data.field(20,20)
tcfield([20,20]')
tcfieldnotrunc([20,20]')
smooth_data.field(1,10)
tcfield([1,10]')
tcfieldnotrunc([1,10]')

%% 2D convolution field (with mask)
Dim = [3,3]; nsubj = 20; resadd = 3; FWHM = 8;
params = ConvFieldParams([FWHM, FWHM], resadd);
mu = zeros(Dim); mu(1,1) = 10; mu(1,3) = 10;
mask = true(Dim); mask(2,1) = 0;
lat_data = wfield(mask, nsubj);
tcfield = convfield_t(lat_data, params);
surf(tcfield.field)
title('2D convolution t field')
xlabel('x'); ylabel('y'); zlabel('t field')

%% 2D convolution field (with mask)
Dim = [10,10]; nsubj = 20; resadd = 3; FWHM = 8;
params = ConvFieldParams([FWHM, FWHM], resadd);
mask = true(Dim); mask(4:7,1:5) = 0;
lat_data = wfield(mask, nsubj);
tcfield = convfield_t(lat_data, params);
surf(tcfield.field)
title('2D convolution t field')
xlabel('x'); ylabel('y'); zlabel('t field')

%% 3D convolution t field
Dim = [11,11,11]; nsubj = 50; resadd = 3; FWHM = 3;
params = ConvFieldParams([FWHM, FWHM, FWHM], resadd);
mask = true(Dim);
lat_data = wfield(mask, nsubj);
tcf = convfield_t( lat_data, params );
surf(tcf.field(:,:,30))
title('3D convolution t field slice')
