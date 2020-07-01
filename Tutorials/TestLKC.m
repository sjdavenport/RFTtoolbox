clear all
close all

% Domain
D     = 3;
T     = 20;
nsubj = 25;
dim   = T * ones( [ 1 D ] );

% Kernel
sigma = 1.1; % smoothing parameter
FWHM  = 3;%sigma2FWHM( sigma );
pad   = ceil(4*sigma);

% Theory LKCs
theoryL = LKC_isogauss_theory( FWHM, dim )

% Construction
resadd  = 1;
enlarge = ceil( resadd / 2 );

% Mask
mask = logical( pad_vals( true( [ dim, 1] ), pad) );


%% LKC computation
%--------------------------------------------------------------------------

% Generate lattice data
lat_data = squeeze( wnfield( mask, nsubj ) );

% logical indicating whether the lat_data gets masked before smoothing or
% not,
lat_masked = false;

field   = Mask( convfield_Field( lat_data, FWHM, 0, resadd, lat_masked, enlarge ) );
dfield  = Mask( convfield_Field( lat_data, FWHM, 1, resadd, lat_masked, enlarge ) );
d2field = Mask( convfield_Field( lat_data, FWHM, 2, resadd, lat_masked, enlarge ) );

% Estimate LKCs
L = LKC_voxmfd_est( field, dfield, d2field, [1 1 1] )

%% Timing the parts if the pipeline depending on inclusion of L1
%--------------------------------------------------------------------------
times_fullL1 = NaN * [ 1 1 1 ];
times_partL1 = NaN * [ 1 1 1 ]; % only stationary integral
times_noL1   = NaN * [ 1 1 1 ];

% Generate convolution fields from lattice data
tic
field   = Mask( convfield_Field( lat_data, FWHM, 0, resadd, lat_masked, enlarge ) );
dfield  = Mask( convfield_Field( lat_data, FWHM, 1, resadd, lat_masked, enlarge ) );
times_partL1(1) = toc;
times_noL1(1)   = times_partL1(1);

tic
d2field = Mask( convfield_Field( lat_data, FWHM, 2, resadd, lat_masked, enlarge ) );
times_fullL1(1) = toc;
times_fullL1(1) = times_fullL1(1) + times_noL1(1);

% Compute timing of Riemannian metric and Christoffel symbols computation
tic
voxmfd = Field2VoxManifold( field, dfield, d2field );
times_fullL1(2) = toc;

% Solve the integrals
tic
LKC_est( voxmfd );
times_fullL1(3) = toc;


% Compute timing of Riemannian metric and Christoffel symbols computation
tic
voxmfd = Field2VoxManifold( field, dfield );
times_partL1(2) = toc;
times_noL1(2)   = times_partL1(2);

% Solve the integrals
tic
LKC_est( voxmfd, [ 1 1 0 ] );
times_partL1(3) = toc;

% Solve the integrals
tic
LKC_est( voxmfd, [ 1 0 0 ] );
times_noL1(3) = toc;

%% Show results
%--------------------------------------------------------------------------
[ times_fullL1;
  times_partL1;
  times_noL1 ]