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
[tcfield, xvals_vecs] = convfield_t(lat_data, FWHM, resadd);
lattice_tfield = convfield_t(lat_data, FWHM, 0);
plot(1:nvox, lattice_tfield, 'o-')
hold on
plot(xvals_vecs{1}, tcfield)
title('1D convolution t field')
legend('Convolution field', 'Lattice Evaluation')
xlabel('voxels')
ylabel('t field')

%% 2D convolution field
Dim = [10,10]; nsubj = 20; resadd = 5; FWHM = 2;
lat_data = normrnd(0,1,[Dim,nsubj]);
tcfield = convfield_t(lat_data, FWHM, resadd);
lattice_tfield = convfield_t(lat_data, FWHM, 0);
zmax = max(tcfield(:)); zmin = min(tcfield(:));
subplot(1,2,1); surf(lattice_tfield); zlim([zmin,zmax])
title('2D lattice t field')
xlabel('x'); ylabel('y'); zlabel('t field')
subplot(1,2,2); surf(tcfield); zlim([zmin,zmax])
title('2D convolution t field')
xlabel('x'); ylabel('y'); zlabel('t field')

%% 2D convolution field (with mask)
Dim = [3,3]; nsubj = 20; resadd = 3; FWHM = 8;
mu = zeros(Dim); mu(1,1) = 10; mu(1,3) = 10;
mask = ones(Dim); mask(2,1) = 0;
lat_data = normrnd(0,1,[Dim,nsubj]).*mask + mu;
mask_hr = zero2nan(mask_highres(logical(mask), resadd, 0) );
tcfield = convfield_t(lat_data, FWHM, resadd);
lattice_tfield = convfield_t(lat_data, FWHM, 0);
surf(tcfield.*mask_hr)
title('2D convolution t field')
xlabel('x'); ylabel('y'); zlabel('t field')


%% 2D convolution field (with mask)
Dim = [10,10]; nsubj = 20; resadd = 3; FWHM = 8;
mask = ones(Dim); mask(4:7,1:5) = 0;
lat_data = normrnd(0,1,[Dim,nsubj]).*mask; 
mask_hr = zero2nan(mask_highres(logical(mask), resadd, 0));
tcfield = convfield_t(lat_data, FWHM, resadd);
lattice_tfield = convfield_t(lat_data, FWHM, 0);
surf(tcfield.*mask_hr)
title('2D convolution t field')
xlabel('x'); ylabel('y'); zlabel('t field')

