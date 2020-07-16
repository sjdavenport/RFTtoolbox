clear all
close all

%% 2D examples
%--------------------------------------------------------------------------
% General parameters
sigma = 5;
T     = 49;
FWHM  = sigma2FWHM( 5 );
nsubj = 10;
pad   = ceil( 4 * sigma );
dim   = [ T T ];
dimp  = dim + 2 * pad;
resadd  = 1;
enlarge = ceil( resadd / 2 );
Mboot   = 500;

%%% Two different voxelmanifolds
% Cube manifold
mask  = true( dim );
mask  = logical( pad_vals( mask, pad) );
cubeL = LKC_isogauss_theory( FWHM, [ T T ] );

% Sphere manifold
mask_sphere = true( dim );
mask_sphere = bndry_voxels( logical( pad_vals( mask_sphere, pad) ), 'full' );

lambda  = 1 / ( 2 * FWHM2sigma( FWHM )^2 );
sphereL = [ sqrt( lambda ), lambda ] .* [ 4*(T-1) 4*(T-1) ];

%% 2D LKC simulation cube
%--------------------------------------------------------------------------
Msim = 50;

L_conv_ests = NaN * ones( [ Msim 2 ] );
L_HP_ests   = NaN * ones( [ Msim 2 ] );
L_bHP_ests  = NaN * ones( [ Msim 2 ] );

tic
for m = 1:Msim
    lat_data = wnfield( mask, nsubj );

    % logical indicating whether the lat_data gets masked before smoothing or
    % not,
    lat_masked = false;

    % Generate convolution fields from lattice data
    cfield  = convfield_Field( lat_data, FWHM, 0, resadd, lat_masked, enlarge );
    dcfield = convfield_Field( lat_data, FWHM, 1, resadd, lat_masked, enlarge );
    
    % Compute the different LKC estimates
    L_conv_ests(m,:) = LKC_voxmfd_est( cfield, dcfield );
    tmp = LKC_HP_est( cfield, 1, 1 );
    L_HP_ests(m,:)   = tmp.hatL;
    tmp = LKC_HP_est( cfield, Mboot, 1 );
    L_bHP_ests(m,:)  = tmp.hatL;
end
sim_time = toc

% Show results
results = ...
    struct( 'theory',    cubeL, ...
            'convE',     mean( L_conv_ests, 1 ), ...
            'std_convE', 1.96*std( L_conv_ests, 0, 1 ) / sqrt( Msim ), ...
            'HPE',       mean( L_HP_ests, 1 ), ...
            'std_HPE',   1.96*std( L_HP_ests, 0, 1 ) / sqrt( Msim ), ...
            'bHPE',      mean( L_bHP_ests, 1 ), ...
            'std_bHPE',  1.96*std( L_bHP_ests, 0, 1 ) / sqrt( Msim ) )
        
% ~ 126.5465 seconds, Msim = 50, nsubj = 10, Mboot = 5e2, resadd = 1,
% FWHM = 11.7741        
%        theory: [13.8593 48.0200]
%         convE: [13.9888 47.8210]
%     std_convE: [0.1913 1.1655]
%           HPE: [13.4792 42.4656]
%       std_HPE: [0.2694 1.1559]
%          bHPE: [14.0071 48.0773]
%      std_bHPE: [0.2185 1.3918]

%% 2D LKC simulation sphere
%--------------------------------------------------------------------------
Msim = 50;

L_conv_ests = NaN * ones( [ Msim 2 ] );
L_HP_ests   = NaN * ones( [ Msim 2 ] );
L_bHP_ests  = NaN * ones( [ Msim 2 ] );

tic
for m = 1:Msim
    lat_data = wnfield( mask_sphere, nsubj );

    % logical indicating whether the lat_data gets masked before smoothing or
    % not,
    lat_masked = false;

    % Generate convolution fields from lattice data
    cfield  = convfield_Field( lat_data, FWHM, 0, resadd, lat_masked, enlarge );
    dcfield = convfield_Field( lat_data, FWHM, 1, resadd, lat_masked, enlarge );

    % Compute the different LKC estimates
    L_conv_ests(m,:) = LKC_voxmfd_est( cfield, dcfield );
    tmp = LKC_HP_est( cfield, 1, 1 );
    L_HP_ests(m,:)   = tmp.hatL;
    tmp = LKC_HP_est( cfield, Mboot, 1 );
    L_bHP_ests(m,:)  = tmp.hatL;
end
sim_time = toc

% Show results
results = ...
    struct( 'theory',    sphereL, ...
            'convE',     mean( L_conv_ests, 1 ), ...
            'std_convE', 1.96*std( L_conv_ests, 0, 1 ) / sqrt( Msim ), ...
            'HPE',       mean( L_HP_ests, 1 ), ...
            'std_HPE',   1.96*std( L_HP_ests, 0, 1 ) / sqrt( Msim ), ...
            'bHPE',      mean( L_bHP_ests, 1 ), ...
            'std_bHPE',  1.96*std( L_bHP_ests, 0, 1 ) / sqrt( Msim ) )
        
% ~ 127.3411 seconds, Msim = 50, nsubj = 10, Mboot = 5e2, resadd = 1,
% FWHM = 11.7741        
%        theory: [27.1529 3.8400]
%         convE: [27.1531 3.9105]
%     std_convE: [0.4225 0.1470]
%           HPE: [26.3976 3.8981]
%       std_HPE: [0.4203 0.5437]
%          bHPE: [27.1322 4.1257]
%      std_bHPE: [0.4447 0.2879]


%% 2D LKC simulation sphere (nonstationary)
%--------------------------------------------------------------------------
Msim = 50;

L_conv_ests = NaN * ones( [ Msim 2 ] );
L_HP_ests   = NaN * ones( [ Msim 2 ] );
L_bHP_ests  = NaN * ones( [ Msim 2 ] );

tic % ~55 seconds
for m = 1:Msim
    lat_data = wnfield( mask_sphere, nsubj );

    % logical indicating whether the lat_data gets masked before smoothing or
    % not,
    lat_masked = true;

    % Generate convolution fields from lattice data
    cfield  = convfield_Field( lat_data, FWHM, 0, resadd, lat_masked, enlarge );
    dcfield = convfield_Field( lat_data, FWHM, 1, resadd, lat_masked, enlarge );

    % Compute the different LKC estimates
    L_conv_ests(m,:) = LKC_voxmfd_est( cfield, dcfield );
    tmp = LKC_HP_est( cfield, 1, 1 );
    L_HP_ests(m,:)   = tmp.hatL;
    tmp = LKC_HP_est( cfield, Mboot, 1 );
    L_bHP_ests(m,:)  = tmp.hatL;
end
sim_time = toc

% Show results
results = ...
    struct( 'theory',    sphereL, ...
            'convE',     mean( L_conv_ests, 1 ), ...
            'std_convE', 1.96*std( L_conv_ests, 0, 1 ) / sqrt( Msim ), ...
            'HPE',       mean( L_HP_ests, 1 ), ...
            'std_HPE',   1.96*std( L_HP_ests, 0, 1 ) / sqrt( Msim ), ...
            'bHPE',      mean( L_bHP_ests, 1 ), ...
            'std_bHPE',  1.96*std( L_bHP_ests, 0, 1 ) / sqrt( Msim ) )
        
% ~132 seconds, Msim = 50, nsubj = 10, Mboot = 5e2, resadd = 1,
% FWHM = 11.7741        
%        theory: [27.1529 3.8400]
%         convE: [26.2490 0.4943]
%     std_convE: [0.4413 0.0303]
%           HPE: [28.1253 0.9700]
%       std_HPE: [0.6902 0.8855]
%          bHPE: [29.0909 0.6848]
%      std_bHPE: [0.5076 0.3093]

% ~266 seconds, Msim = 50, nsubj = 10, Mboot = 5e2, resadd = 1,
% FWHM = 11.7741
%        theory: [27.1529 3.8400]
%         convE: [26.2602 0.2515]
%     std_convE: [0.3228 0.0158]
%           HPE: [28.3852 -0.4676]
%       std_HPE: [0.5018 0.8810]
%          bHPE: [29.1521 0.2722]
%      std_bHPE: [0.3507 0.3519]


%% 3D example
%--------------------------------------------------------------------------
clear all
close all

% General parameters
D     = 3;
T     = 5;
FWHM  = 2;
nsubj = 10;
pad   = ceil( 4 * FWHM2sigma( FWHM ) );
dim   = T * ones( [ 1 D ] );
dimp  = dim + 2 * pad;
resadd  = 3;
enlarge = ceil( resadd / 2 );
xvals   = { 1:dimp(1), 1:dimp(2) };
Mboot   = 100;

%%% Two different voxelmanifolds
% Cube manifold
mask = true( dim );
mask = logical( pad_vals( mask, pad ) );
cubeL = LKC_isogauss_theory( FWHM, [ T T T ] );

% Sphere manifold
mask_sphere = true( dim );
mask_sphere = bndry_voxels( logical( pad_vals( mask_sphere, pad ) ), 'full' );
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

L_conv_ests = NaN * ones( [ Msim 3 ] );
L_HP_ests   = NaN * ones( [ Msim 3 ] );
L_bHP_ests  = NaN * ones( [ Msim 3 ] );

tic % ~60 seconds
for m = 1:Msim
%     m
    lat_data = wnfield( mask, nsubj );

    % logical indicating whether the lat_data gets masked before smoothing or
    % not,
    lat_masked = false;

    % Generate convolution fields from lattice data
    cfield   = convfield_Field( lat_data, FWHM, 0, resadd, lat_masked, enlarge );
    dcfield  = convfield_Field( lat_data, FWHM, 1, resadd, lat_masked, enlarge );
    d2cfield = convfield_Field( lat_data, FWHM, 2, resadd, lat_masked, enlarge );
    
    L_conv_ests(m,:) = LKC_voxmfd_est( cfield, dcfield, d2cfield );
    tmp = LKC_HP_est( cfield, 1, 1 );
    L_HP_ests(m,:)   = tmp.hatL;
    tmp = LKC_HP_est( cfield, Mboot, 1 );
    L_bHP_ests(m,:)  = tmp.hatL;
end
toc

% Show results
results = ...
    struct( 'theory',    cubeL, ...
            'convE',     mean( L_conv_ests, 1 ), ...
            'std_convE', 1.96*std( L_conv_ests, 0, 1 ) / sqrt( Msim ), ...
            'HPE',       mean( L_HP_ests, 1 ), ...
            'std_HPE',   1.96*std( L_HP_ests, 0, 1 ) / sqrt( Msim ), ...
            'bHPE',      mean( L_bHP_ests, 1 ), ...
            'std_bHPE',  1.96*std( L_bHP_ests, 0, 1 ) / sqrt( Msim ) )
        
% ~ 765 seconds, Msim = 50, nsubj = 10, Mboot = 2e2, resadd = 3,
% FWHM = 3
%        theory: [8.3255 23.1049 21.3734]
%         convE: [8.0028 23.1785 21.0434]
%     std_convE: [0.3007 0.5363 1.2732]
%           HPE: [8.1326 20.5514 11.7887]
%       std_HPE: [0.1929 0.6420 1.2553]
%          bHPE: [8.3129 23.0012 20.2729]
%      std_bHPE: [0.1598 0.6057 1.4727]

%% 3D LKC estimation cube manifold L1 only second integral (nonstationary part)
%--------------------------------------------------------------------------

Msim = 50;

L_conv_ests = NaN * ones( [ Msim 3 ] );
L_HP_ests   = NaN * ones( [ Msim 3 ] );
L_bHP_ests  = NaN * ones( [ Msim 3 ] );

tic
for m = 1:Msim
%     m
    lat_data = wnfield( mask, nsubj );

    % logical indicating whether the lat_data gets masked before smoothing or
    % not,
    lat_masked = false;

    % Generate convolution fields from lattice data
    cfield   = convfield_Field( lat_data, FWHM, 0, resadd, lat_masked, enlarge );
    dcfield  = convfield_Field( lat_data, FWHM, 1, resadd, lat_masked, enlarge );
    d2cfield = convfield_Field( lat_data, FWHM, 2, resadd, lat_masked, enlarge );
    
    L_conv_ests(m,:) = LKC_voxmfd_est( cfield, dcfield, d2cfield, [1 0 1] );
    tmp = LKC_HP_est( cfield, 1, 1 );
    L_HP_ests(m,:)   = tmp.hatL;
    tmp = LKC_HP_est( cfield, Mboot, 1 );
    L_bHP_ests(m,:)  = tmp.hatL;
end
toc

% Show results
results = ...
    struct( 'theory', cubeL, ...
            'convE', mean( L_conv_ests, 1 ), ...
            'std_convE', 1.96*std( L_conv_ests, 0, 1 ) / sqrt( Msim ), ...
            'HPE', mean( L_HP_ests, 1 ), ...
            'std_HPE', 1.96*std( L_HP_ests, 0, 1 ) / sqrt( Msim ), ...
            'bHPE', mean( L_bHP_ests, 1 ), ...
            'std_bHPE', 1.96*std( L_bHP_ests, 0, 1 ) / sqrt( Msim ) )

%% 3D LKC estimation sphere manifold all integrals (stationary)
%--------------------------------------------------------------------------

Msim = 20;

L_conv_ests = NaN * ones( [ Msim 3 ] );
L_HP_ests   = NaN * ones( [ Msim 3 ] );
L_bHP_ests  = NaN * ones( [ Msim 3 ] );

tic
for m = 1:Msim
%     m
    lat_data = wnfield( mask_sphere, nsubj );

    % logical indicating whether the lat_data gets masked before smoothing or
    % not,
    lat_masked = false;

    % Generate convolution fields from lattice data
    cfield   = convfield_Field( lat_data, FWHM, 0, resadd, lat_masked, enlarge );
    dcfield  = convfield_Field( lat_data, FWHM, 1, resadd, lat_masked, enlarge );
    d2cfield = convfield_Field( lat_data, FWHM, 2, resadd, lat_masked, enlarge );
    
    L_conv_ests(m,:) = LKC_voxmfd_est( cfield, dcfield, d2cfield, [1 1 1] );
    tmp = LKC_HP_est( cfield, 1, 1 );
    L_HP_ests(m,:)   = tmp.hatL;
    tmp = LKC_HP_est( cfield, Mboot, 1 );
    L_bHP_ests(m,:)  = tmp.hatL;
end
toc

% Show results
results = ...
    struct( 'theory',    sphereL, ...
            'convE',     mean( L_conv_ests, 1 ), ...
            'std_convE', 1.96*std( L_conv_ests, 0, 1 ) / sqrt( Msim ), ...
            'HPE',       mean( L_HP_ests, 1 ), ...
            'std_HPE',   1.96*std( L_HP_ests, 0, 1 ) / sqrt( Msim ), ...
            'bHPE',      mean( L_bHP_ests, 1 ), ...
            'std_bHPE',  1.96*std( L_bHP_ests, 0, 1 ) / sqrt( Msim ) )

% ~ 95 seconds, Msim = 20, nsubj = 10, Mboot = 2e2, resadd = 1,
% FWHM = 5 
%        theory: [1.9981 11.3122 3.6195]
%         convE: [2.0504 11.6608 3.6769]
%     std_convE: [0.2799 1.0503 0.5425]
%           HPE: [2.0114 9.5706 1.6963]
%       std_HPE: [0.1695 0.9739 0.9715]
%          bHPE: [2.1421 11.7316 3.8366]
%      std_bHPE: [0.1589 1.2156 1.0306]   

% ~ 50 seconds, Msim = 20, nsubj = 10, Mboot = 2e2, resadd = 1,
% FWHM = 3 
%        theory: [3.3302 31.4227 16.7568]
%         convE: [3.8013 31.8278 17.6865]
%     std_convE: [0.2841 1.2554 1.2824]
%           HPE: [3.1921 27.2940 10.7415]
%       std_HPE: [0.3007 1.1808 1.6568]
%          bHPE: [3.3222 30.7598 15.9880]
%      std_bHPE: [0.2495 1.5999 1.9516]

% ~ 305 seconds, Msim = 20, nsubj = 10, Mboot = 2e2, resadd = 3,
% FWHM = 3 
%        theory: [3.3302 31.4227 16.7568]
%         convE: [3.6070 31.5453 16.3614]
%     std_convE: [0.4178 1.3974 1.2814]
%           HPE: [3.2962 27.9853 10.5255]
%       std_HPE: [0.4489 1.5299 1.4893]
%          bHPE: [3.3569 31.5809 16.2687]
%      std_bHPE: [0.2770 1.7259 1.9246]