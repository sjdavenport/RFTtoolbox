clear all
close all

%% 2D examples
%--------------------------------------------------------------------------
% General parameters
sigma = 5;
T     = 49;
FWHM  = sigma2FWHM( 5 );
nsubj = 10;
pad   = ceil(4*sigma);
dim   = [ T T ];
dimp  = dim + 2 * pad;
resadd  = 1;
enlarge = ceil( resadd / 2 );
xvals   = { (1:dimp(1)), (1:dimp(2)) };

%%% Two different voxelmanifolds
% Cube manifold
mask  = true( dim );
mask  = logical( pad_vals( mask, pad) );
cubeL = LKC_isogauss_theory( FWHM, [ T T ] )

% Sphere manifold
mask_sphere = true( dim );
mask_sphere = bndry_voxels( logical( pad_vals( mask_sphere, pad) ), 'full' );

lambda  = 1 / ( 2 * FWHM2sigma( FWHM )^2 );
sphereL = [ sqrt( lambda ), lambda ] .* [ 4*(T-1) 4*(T-1) ]

%% 2D LKC simulation cube
%--------------------------------------------------------------------------
Msim = 50;

Lests = NaN * ones( [ Msim 2 ] );

tic % couple of seconds
for m = 1:Msim
    lat_data = wnfield( mask, nsubj );

    % logical indicating whether the lat_data gets masked before smoothing or
    % not,
    lat_masked = false;

    % Generate convolution fields from lattice data
    cfield  = convfield_Field( lat_data, FWHM, 0, resadd, lat_masked, enlarge );
    dcfield = convfield_Field( lat_data, FWHM, 1, resadd, lat_masked, enlarge );
    
    Lests(m,:) = LKC_voxmfd_est( cfield, dcfield );
end
toc

% 
[ mean(Lests, 1);
  1.96*std( Lests, 0, 1 ) / sqrt( Msim );
  cubeL ]

%% 2D LKC simulation sphere
%--------------------------------------------------------------------------
Msim = 50;

Lests = NaN * ones( [ Msim 2 ] );

tic % couple of seconds
for m = 1:Msim
    lat_data = wnfield( mask_sphere, nsubj );

    % logical indicating whether the lat_data gets masked before smoothing or
    % not,
    lat_masked = false;

    % Generate convolution fields from lattice data
    cfield  = convfield_Field( lat_data, FWHM, 0, resadd, lat_masked, enlarge );
    dcfield = convfield_Field( lat_data, FWHM, 1, resadd, lat_masked, enlarge );
    
    Lests(m,:) = LKC_voxmfd_est( cfield, dcfield );
end
toc

% 
[ mean(Lests, 1);
  1.96*std( Lests, 0, 1 ) / sqrt( Msim );
  sphereL ]

%% 3D example
%--------------------------------------------------------------------------
clear all
close all
% General parameters
D     = 3;
T     = 5;
FWHM  = 5;
nsubj = 10;
pad   = ceil( 4 * FWHM2sigma( FWHM ) );
dim   = T * ones( [ 1 D ] );
dimp  = dim + 2 * pad;
resadd  = 3;
enlarge = ceil( resadd / 2 );
xvals   = { (1:dimp(1)), (1:dimp(2)) };

%%% Two different voxelmanifolds
% Cube manifold
mask = true( dim );
mask = logical( pad_vals( mask, pad) );
cubeL = LKC_isogauss_theory( FWHM, [ T T T ] );

% Sphere manifold
mask_sphere = true( dim );
mask_sphere = bndry_voxels( logical( pad_vals( mask_sphere, pad) ), 'full' );
sm = size(mask_sphere);

lambda   = 1 / ( 2 * FWHM2sigma( FWHM )^2 );
volmask  = 2*T^2 + 2*(T-2)*T + 2* (T-2)^2;
surfmask = 6 * ( T^2 + (T-2)^2 );
linemask = 12 * T / 4 - 12 * ( T - 2 ) * 1 / 4;

sphere = VoxManifold( mask_highres( mask_sphere, 1 ), ...
                      xvals_highres( {1:sm(1),1:sm(2),1:sm(3)}, 1 ) );
sphereL  = [ sqrt(lambda), lambda, lambda^(3/2) ] .* LKC_est( sphere );
sphereL1 = [ sqrt(lambda), lambda, lambda^(3/2) ] .* [ NaN surfmask/2 volmask ];

%% 3D LKC estimation cube manifold all integrals
%--------------------------------------------------------------------------
Msim = 50;

Lests = NaN * ones( [ Msim D ] );
tic % ~60 seconds
for m = 1:Msim
    m
    lat_data = wnfield( mask, nsubj );

    % logical indicating whether the lat_data gets masked before smoothing or
    % not,
    lat_masked = false;

    % Generate convolution fields from lattice data
    cfield   = convfield_Field( lat_data, FWHM, 0, resadd, lat_masked, enlarge );
    dcfield  = convfield_Field( lat_data, FWHM, 1, resadd, lat_masked, enlarge );
    d2cfield = convfield_Field( lat_data, FWHM, 2, resadd, lat_masked, enlarge );
    
    Lests(m,:) = LKC_voxmfd_est( cfield, dcfield, d2cfield );
end
toc

% 
[ mean(Lests, 1);...
  cubeL;...
  1.96*std( Lests, 0, 1 )/sqrt(Msim)
  ]

%% 3D LKC estimation cube manifold L1 only second integral (nonstationary part)
%--------------------------------------------------------------------------
Msim = 50;

Lests = NaN * ones( [ Msim D ] );
tic % ~2min
for m = 1:Msim
    m
    lat_data = wnfield( mask, nsubj );

    % logical indicating whether the lat_data gets masked before smoothing or
    % not,
    lat_masked = false;

    % Generate convolution fields from lattice data
    cfield   = convfield_Field( lat_data, FWHM, 0, resadd, lat_masked, enlarge );
    dcfield  = convfield_Field( lat_data, FWHM, 1, resadd, lat_masked, enlarge );
    d2cfield = convfield_Field( lat_data, FWHM, 2, resadd, lat_masked, enlarge );
    
    Lests(m,:) = LKC_voxmfd_est( cfield, dcfield, d2cfield, [1 0 1] );
end
toc

% Note the main part in L1 is ignored. Hence L1 should be close to zero
[ mean( Lests, 1 );...
  cubeL;...
  1.96 * std( Lests, 0, 1 ) / sqrt( Msim )
 ]

%% 3D LKC estimation sphere manifold all integrals
%--------------------------------------------------------------------------
Msim = 20;

Lests = NaN * ones( [ Msim D ] );
tic % ~60 seconds
for m = 1:Msim
    m
    lat_data = wnfield( mask_sphere, nsubj );

    % logical indicating whether the lat_data gets masked before smoothing or
    % not,
    lat_masked = false;

    % Generate convolution fields from lattice data
    cfield   = convfield_Field( lat_data, FWHM, 0, resadd, lat_masked, enlarge );
    dcfield  = convfield_Field( lat_data, FWHM, 1, resadd, lat_masked, enlarge );
    d2cfield = convfield_Field( lat_data, FWHM, 2, resadd, lat_masked, enlarge );
    
    Lests(m,:) = LKC_voxmfd_est( cfield, dcfield, d2cfield, [1 1 1] );
end
toc

% 
[ mean(Lests, 1);...
  sphereL;...
  1.96 * std( Lests, 0, 1 ) / sqrt( Msim )
  ]