clear all
close all

%% 2D example
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
theoryL = LKC_isogauss_theory( FWHM, [ T T ] )

% Mask
mask = true( dim );
mask = logical( pad_vals( mask, pad) );

%% Generate lattice data and get convolution fields small simulation
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
  theoryL ]


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
resadd  = 1;
enlarge = ceil( resadd / 2 );
xvals   = { (1:dimp(1)), (1:dimp(2)) };
theoryL = LKC_isogauss_theory( FWHM, [ T T T ] )

% Mask
mask = true( dim );
mask = logical( pad_vals( mask, pad) );

%% Generate lattice data and get convolution fields small simulation for
%  ignoring nonstationary integral in LKC estimation 
%--------------------------------------------------------------------------
Msim = 50;

Lests = NaN * ones( [ Msim D ] );
tic % ~60 seconds
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
[ mean(Lests, 1);...
  theoryL;...
  1.96*std( Lests, 0, 1 )/sqrt(Msim)
  ]

%% Generate lattice data and get convolution fields small simulation for
%  ignoring nonstationary integral in LKC estimation 
%--------------------------------------------------------------------------
Msim = 50;

Lests = NaN * ones( [ Msim D ] );
tic % ~2min
for m = 1:Msim
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
  theoryL;...
  1.96*std( Lests, 0, 1 )/sqrt(Msim)
 ]