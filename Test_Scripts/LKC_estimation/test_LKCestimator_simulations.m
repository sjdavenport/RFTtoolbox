%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the different LKC estimates with short simualtions
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D examples
%--------------------------------------------------------------------------
% General parameters
D     = 1; 
sigma = 2;
T     = 100;
FWHM  = sigma2FWHM( sigma );
nsubj = 10;
pad   = ceil( 4 * sigma );
dim   = [ T 1 ];
dimp  = dim + 2 * [ pad, 0 ];
resadd  = 1;
enlarge = ceil( resadd / 2 );
Mboot   = 200;
Msim    = 50;

%%% Two different voxelmanifolds
% Cube manifold
mask  = true( dim );
mask  = logical( pad_vals( mask, pad) );
cubeL = LKC_isogauss_theory( FWHM, dim );

% Sphere manifold
mask2 = true( dim );
mask2( 30:60 ) = false;
mask2 = logical( pad_vals( mask2, pad) );

lambda  = 1 / ( 2 * FWHM2sigma( FWHM )^2 );
Lmask2  = sqrt(lambda) * sum(mask2);

%% 1D LKC simulation cube
%--------------------------------------------------------------------------

params = ConvFieldParams( FWHM*ones( [ 1 D ] ), resadd, enlarge, false );
methods = struct( 'convE', true, 'HPE', true, 'bHPE', [ Mboot, 1 ] );

results = simulate_LKCests( Msim, nsubj, methods, params,...
                                     mask );
% Show results
struct( 'theory',    cubeL, ...
        'convE',     mean( results.convE, 1 ), ...
        'std_convE', 1.96*std( results.convE, 0, 1 ) / sqrt( Msim ), ...
        'HPE',       mean( results.HPE, 1 ), ...
        'std_HPE',   1.96*std( results.HPE, 0, 1 ) / sqrt( Msim ), ...
        'bHPE',      mean( results.bHPE, 1 ), ...
        'std_bHPE',  1.96*std( results.bHPE, 0, 1 ) / sqrt( Msim ) )

%% 1D LKC simulation mask2
%--------------------------------------------------------------------------

params = ConvFieldParams( FWHM*ones( [ 1 D ] ), resadd, enlarge, false );
methods = struct( 'convE', true, 'HPE', true, 'bHPE', [ Mboot, 1 ] );

results = simulate_LKCests( Msim, nsubj, methods, params,...
                                     mask2 );
% Show results
struct( 'theory',    Lmask2, ...
        'convE',     mean( results.convE, 1 ), ...
        'std_convE', 1.96*std( results.convE, 0, 1 ) / sqrt( Msim ), ...
        'HPE',       mean( results.HPE, 1 ), ...
        'std_HPE',   1.96*std( results.HPE, 0, 1 ) / sqrt( Msim ), ...
        'bHPE',      mean( results.bHPE, 1 ), ...
        'std_bHPE',  1.96*std( results.bHPE, 0, 1 ) / sqrt( Msim ) )
    

%% %% 2D examples
%--------------------------------------------------------------------------
clear all
close all

% General parameters
D     = 2;
sigma = 5;
T     = 49;
FWHM  = sigma2FWHM( 5 );
nsubj = 10;
pad   = ceil( 4 * sigma );
dim   = T * ones( [ 1 D ] );
dimp  = dim + 2 * pad;
resadd  = 1;
enlarge = ceil( resadd / 2 );
Mboot   = 200;
Msim    = 20;

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

params = ConvFieldParams( FWHM*ones( [ 1 D ] ), resadd, enlarge,...
                          false );
methods = struct( 'convE', true, 'HPE', true, 'bHPE', [ Mboot, 1 ] );

results = simulate_LKCests( Msim, nsubj, methods, params,...
                                     mask );
% Show results
struct( 'theory',    cubeL, ...
        'convE',     mean( results.convE, 1 ), ...
        'std_convE', 1.96*std( results.convE, 0, 1 ) / sqrt( Msim ), ...
        'HPE',       mean( results.HPE, 1 ), ...
        'std_HPE',   1.96*std( results.HPE, 0, 1 ) / sqrt( Msim ), ...
        'bHPE',      mean( results.bHPE, 1 ), ...
        'std_bHPE',  1.96*std( results.bHPE, 0, 1 ) / sqrt( Msim ) )

    
%% 2D LKC simulation sphere
%--------------------------------------------------------------------------

params = ConvFieldParams( FWHM*ones( [ 1 D ] ), resadd, enlarge,...
                          false );
methods = struct( 'convE', true, 'HPE', true, 'bHPE', [ Mboot, 1 ] );

results = simulate_LKCests( Msim, nsubj, methods, params,...
                            mask_sphere );
% Show results
struct( 'theory',    sphereL, ...
        'convE',     mean( results.convE, 1 ), ...
        'std_convE', 1.96*std( results.convE, 0, 1 ) / sqrt( Msim ), ...
        'HPE',       mean( results.HPE, 1 ), ...
        'std_HPE',   1.96*std( results.HPE, 0, 1 ) / sqrt( Msim ), ...
        'bHPE',      mean( results.bHPE, 1 ), ...
        'std_bHPE',  1.96*std( results.bHPE, 0, 1 ) / sqrt( Msim ) )
        

%% 2D LKC simulation sphere (nonstationary)
%--------------------------------------------------------------------------

params = ConvFieldParams( FWHM*ones( [ 1 D ] ), resadd, enlarge,...
                          true );
methods = struct( 'convE', true, 'HPE', true, 'bHPE', [ Mboot, 1 ] );

results = simulate_LKCests( Msim, nsubj, methods, params,...
                                     mask );
% Show results
struct( 'theory',    sphereL, ...
        'convE',     mean( results.convE, 1 ), ...
        'std_convE', 1.96*std( results.convE, 0, 1 ) / sqrt( Msim ), ...
        'HPE',       mean( results.HPE, 1 ), ...
        'std_HPE',   1.96*std( results.HPE, 0, 1 ) / sqrt( Msim ), ...
        'bHPE',      mean( results.bHPE, 1 ), ...
        'std_bHPE',  1.96*std( results.bHPE, 0, 1 ) / sqrt( Msim ) )

        
%% %% 3D examples
%--------------------------------------------------------------------------
% Clean workspace
clear all
close all

% General parameters
D     = 3;
T     = 5;
FWHM  = 4;
nsubj = 10;
pad   = ceil( 4 * FWHM2sigma( FWHM ) );
dim   = T * ones( [ 1 D ] );
dimp  = dim + 2 * pad;
resadd  = 1;
enlarge = ceil( resadd / 2 );
xvals   = { 1:dimp(1), 1:dimp(2), 1:dimp(3) };
Mboot   = 100;
Msim    = 20;

%%% Two different voxelmanifolds
% Cube manifold
mask  = true( dim );
mask  = logical( pad_vals( mask, pad ) );
cubeL = LKC_isogauss_theory( FWHM, [ T T T ] );

% Sphere manifold
mask_sphere = true( dim );
mask_sphere = bndry_voxels( logical( pad_vals( mask_sphere, pad ) ), 'full' );
sm = size( mask_sphere );

lambda   = 1 / ( 2 * FWHM2sigma( FWHM )^2 );
volmask  = 2*T^2 + 2*(T-2)*T + 2* (T-2)^2;
surfmask = 6 * ( T^2 + (T-2)^2 );
linemask = 12 * T / 4 - 12 * ( T - 2 ) * 1 / 4;

sphere = VoxManifold( mask_highres( mask_sphere, 1 ), ...
                      xvals_highres( {1:sm(1),1:sm(2),1:sm(3)}, 1 ) );
sphereL  = [ sqrt(lambda), lambda, lambda^(3/2) ] .* LKC_est( sphere );
sphereL1 = [ sqrt(lambda), lambda, lambda^(3/2) ] .* [ NaN surfmask/2 volmask ];

%% 3D cube manifold (stationary)
%--------------------------------------------------------------------------

params = ConvFieldParams( FWHM*ones( [ 1 D ] ), resadd, enlarge,...
                          false );
methods = struct( 'convE', [ true, true, true ], 'HPE', true, 'bHPE', [ Mboot, 1 ] );

results = simulate_LKCests( Msim, nsubj, methods, params, mask );

% Show results
struct( 'theory',    cubeL, ...
        'convE',     mean( results.convE, 1 ), ...
        'std_convE', 1.96*std( results.convE, 0, 1 ) / sqrt( Msim ), ...
        'L1_locstat', mean( results.convE(:,1)-results.nonstatInt ), ...
        'std_L1_locstat', std( results.convE(:,1)-results.nonstatInt ), ...
        'HPE',       mean( results.HPE, 1 ), ...
        'std_HPE',   1.96*std( results.HPE, 0, 1 ) / sqrt( Msim ), ...
        'bHPE',      mean( results.bHPE, 1 ), ...
        'std_bHPE',  1.96*std( results.bHPE, 0, 1 ) / sqrt( Msim ) )
        
%% 3D sphere manifold (stationary)
%--------------------------------------------------------------------------

params = ConvFieldParams( FWHM*ones( [ 1 D ] ), resadd, enlarge,...
                          false );
methods = struct( 'convE', [ true, true, true ], 'HPE', true, 'bHPE', [ Mboot, 1 ] );

results = simulate_LKCests( Msim, nsubj, methods, params, mask_sphere );

% Show results
struct( 'theory',    cubeL, ...
        'convE',     mean( results.convE, 1 ), ...
        'std_convE', 1.96*std( results.convE, 0, 1 ) / sqrt( Msim ), ...
        'L1_locstat', mean( results.convE(:,1)-results.nonstatInt ), ...
        'std_L1_locstat', 1.96*std( results.convE(:,1)-results.nonstatInt ) / sqrt( Msim ), ...
        'HPE',       mean( results.HPE, 1 ), ...
        'std_HPE',   1.96*std( results.HPE, 0, 1 ) / sqrt( Msim ), ...
        'bHPE',      mean( results.bHPE, 1 ), ...
        'std_bHPE',  1.96*std( results.bHPE, 0, 1 ) / sqrt( Msim ) )

%% 3D sphere manifold (nonstationary)
%--------------------------------------------------------------------------

params = ConvFieldParams( FWHM*ones( [ 1 D ] ), resadd, enlarge,...
                          true );
methods = struct( 'convE', [ true, true, true ], 'HPE', true, 'bHPE', [ Mboot, 1 ] );

results = simulate_LKCests( Msim, nsubj, methods, params, mask_sphere );

% Show results
struct( 'theory',    cubeL, ...
        'convE',     mean( results.convE, 1 ), ...
        'std_convE', 1.96*std( results.convE, 0, 1 ) / sqrt( Msim ), ...
        'L1_locstat', mean( results.convE(:,1)-results.nonstatInt ), ...
        'std_L1_locstat', 1.96*std( results.convE(:,1)-results.nonstatInt ) / sqrt( Msim ), ...
        'HPE',       mean( results.HPE, 1 ), ...
        'std_HPE',   1.96*std( results.HPE, 0, 1 ) / sqrt( Msim ), ...
        'bHPE',      mean( results.bHPE, 1 ), ...
        'std_bHPE',  1.96*std( results.bHPE, 0, 1 ) / sqrt( Msim ) )
    
% Euclidean angles:    
%             theory: [6.2442 12.9965 9.0169]
%              convE: [0.9450 11.6917 2.7323]
%          std_convE: [0.1576 0.7978 0.3643]
%         L1_locstat: 2.0196
%     std_L1_locstat: 0.8587
%                HPE: [1.5129 9.8207 0.9537]
%            std_HPE: [0.1304 0.7102 0.9554]
%               bHPE: [1.4751 11.8202 2.7163]
%           std_bHPE: [0.1576 0.9450 0.9321]

% induced metric angles:    
%             theory: [6.2442 12.9965 9.0169]
%              convE: [1.0458 11.7810 2.7530]
%          std_convE: [0.1022 0.7738 0.2613]
%         L1_locstat: 1.9251
%     std_L1_locstat: 0.3488
%                HPE: [1.3595 10.8083 1.6856]
%            std_HPE: [0.1162 0.8030 0.9676]
%               bHPE: [1.5806 11.6767 3.6200]
%           std_bHPE: [0.1300 0.6733 1.0845]