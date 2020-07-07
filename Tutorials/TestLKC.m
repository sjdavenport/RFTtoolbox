clear all
close all

% Domain
D     = 3;
T     = 5;
nsubj = 25;
dim   = T * ones( [ 1 D ] );

% Kernel
FWHM1  = [ 5 5 5 ];%sigma2FWHM( sigma );
FWHM2  = [ 2 1 3 ];
pad   = ceil( 4 * FWHM2sigma( FWHM1 ) );

% Theory LKCs
theoryL = LKC_isogauss_theory( FWHM1(1), dim )

% Construction
resadd  = 1;
enlarge = ceil( resadd / 2 );

% Mask
mask = logical( pad_vals( true( [ dim, 1] ), pad) );




%% LKC computation
%--------------------------------------------------------------------------

% Generate lattice data
lat_data = squeeze( wnfield( mask, nsubj ) );
lat_data = convfield_Field( lat_data, FWHM1, 0, resadd, false, enlarge );
%M = lat_data.masksize;
%lat_data.field = lat_data.field( randsample( M(1), M(1) )', :, randsample( M(3), M(3) )', : );

% logical indicating whether the lat_data gets masked before smoothing or
% not,
lat_masked = false;

field   = Mask( convfield_Field( lat_data, FWHM1, 0, resadd, lat_masked, enlarge ) );
dfield  = Mask( convfield_Field( lat_data, FWHM1, 1, resadd, lat_masked, enlarge ) );
%d2field = Mask( convfield_Field( lat_data, FWHM1, 2, resadd, lat_masked, enlarge ) );
L1 = LKC_voxmfd_est( field, dfield )
theoryL

%%
% Stationary estimator
voxmfd = Field2VoxManifold( field, dfield, d2field );
LKC_est_stationary( voxmfd )

%lat_data.field = permute( lat_data.field, [2 1 3 4] );
field2   = Mask( convfield_Field( lat_data, FWHM2, 0, resadd, lat_masked, enlarge ) );
dfield2  = Mask( convfield_Field( lat_data, FWHM2, 1, resadd, lat_masked, enlarge ) );
d2field2 = Mask( convfield_Field( lat_data, FWHM2, 2, resadd, lat_masked, enlarge ) );

L1 = LKC_voxmfd_est( field, dfield, d2field, [ 1 1 1 ] )
theoryL
Lint1 = LKC_voxmfd_est( field, dfield, d2field, [ 1 0 1 ] )


%%
L2 = LKC_voxmfd_est( field2, dfield2, d2field2, [ 1 1 1 ] )

L3 = LKC_voxmfd_est( field, dfield, d2field, [ 1 0 1 ] )
L4 = LKC_voxmfd_est( field2, dfield2, d2field2, [ 1 1 0 ] )

%% EEC curve and threshold
EEC1 = EEC( -6:0.1:6, L, L0 );
EEC2 = EEC( -6:0.1:6, [ 0, L(1:2) ], L0 );
EEC3 = EEC( -6:0.1:6, L2, L0 );
EEC4 = EEC( -6:0.1:6, L3, L0 );

figure, clf, hold on
plot( EEC1 )
plot( EEC2 )
plot( EEC3 )
plot( EEC4 )

u1 = EECthreshold( 0.05, L, L0 );
u2 = EECthreshold( 0.05, [ 0, L(1:2) ], L0 );
u3 = EECthreshold( 0.05, L2, L0 );
u4 = EECthreshold( 0.05, L3, L0 );

[ u1, u2, u3, u4 ]

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