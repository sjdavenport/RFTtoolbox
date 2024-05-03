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
T     = 50;
FWHM  = sigma2FWHM( 5 );
nsubj = 10;
pad   = ceil( 4 * sigma );
dim   = T * ones( [ 1 D ] );
dimp  = dim + 2 * pad;
resadd  = 0;
enlarge = ceil( resadd / 2 );
Mboot   = 500;
Msim    = 500;

%%% Two different voxelmanifolds
% Cube manifold
mask  = true( dim );
mask  = logical( pad_vals( mask, pad) );
cubeL = LKC_isogauss_theory( FWHM, [ T-1 T-1 ] );

% Sphere manifold
mask_sphere = true( dim );
mask_sphere = bndry_voxels( logical( pad_vals( mask_sphere, pad) ), 'full' );

lambda  = 1 / ( 2 * FWHM2sigma( FWHM )^2 );
sphereL = [ sqrt( lambda ), lambda ] .* [ 4*(T-1) 4*(T-1) ];


%% 2D LKC simulation cube
%--------------------------------------------------------------------------

params = ConvFieldParams( FWHM*ones( [ 1 D ] ), resadd, enlarge,...
                          false );
methods = struct( 'convE', true, 'HPE', true, 'bHPE', [ Mboot, 1 ], 'warpE', 1 );
methods = struct( 'convE', true, 'HPE', false, 'warpE', 0 );
methods = struct( 'HPE', true, 'bHPE', [ Mboot, 1 ], 'warpE', 1 );


results = simulate_LKCests( Msim, nsubj, methods, params,...
                                     mask );
% Show results
struct( 'theory',    cubeL, ...
        ...%'convE',     mean( results.convE, 1 ), ...
        ...%'std_convE', 1.96*std( results.convE, 0, 1 ) / sqrt( Msim ), ...
        'HPE',       mean( results.HPE, 1 ), ...
        'std_HPE',   1.96*std( results.HPE, 0, 1 ) / sqrt( Msim ), ...
        'bHPE',      mean( results.bHPE, 1 ), ...
        'std_bHPE',  1.96*std( results.bHPE, 0, 1 ) / sqrt( Msim ), ...
        'warpE',      mean( results.warpE, 1 ), ...
        'std_warpE',  1.96*std( results.warpE, 0, 1 ) / sqrt( Msim ) )
  
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
methods = struct( 'convE', [ true, true, false ], 'HPE', true, 'bHPE', [ Mboot, 1 ] );

results = simulate_LKCests( Msim, nsubj, methods, params, mask_sphere );

% Show results
struct( 'theory',    sphereL, ...
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
% Note that this gives good agreement between bHPE and convE, if 'convE',
% [ true, true, false ], if 'convE', [ true, true, true ] there is no good
% agreement suggesting that there is a bug in L1 second integral.
methods = struct( 'convE', [ true, true, false ], 'HPE', true, 'bHPE', [ Mboot, 1 ] );

results = simulate_LKCests( Msim, nsubj, methods, params, mask_sphere );

% Show results
struct( 'theory',    sphereL, ...
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

% current implementation: resadd = 1
%             theory: [2.4977 17.6753 7.0693]
%              convE: [1.5166 12.8532 3.2137]
%          std_convE: [0.1708 0.9250 0.3408]
%         L1_locstat: NaN
%     std_L1_locstat: NaN
%                HPE: [1.4789 10.6890 2.3206]
%            std_HPE: [0.1560 0.8508 1.0461]
%               bHPE: [1.5168 12.0942 3.1890]
%           std_bHPE: [0.1832 1.0002 0.9589]

% current implementation: resadd = 1
%             theory: [2.4977 17.6753 7.0693]
%              convE: [1.5166 12.8532 3.2137]
%          std_convE: [0.1708 0.9250 0.3408]
%         L1_locstat: NaN
%     std_L1_locstat: NaN
%                HPE: [1.4789 10.6890 2.3206]
%            std_HPE: [0.1560 0.8508 1.0461]
%               bHPE: [1.5168 12.0942 3.1890]
%           std_bHPE: [0.1832 1.0002 0.9589]

% resadd = 1
% nsubj = 100
% 
% nsubj =
% 
%    100
% 
% 
% ans = 
% 
%   struct with fields:
% 
%             theory: [2.4977 17.6753 7.0693]
%              convE: [1.4682 11.8680 2.8374]
%          std_convE: [0.0247 0.2284 0.0750]
%         L1_locstat: NaN
%     std_L1_locstat: NaN
%                HPE: [1.4626 11.6357 2.3297]
%            std_HPE: [0.0357 0.2178 0.4184]
%               bHPE: [1.5370 11.5182 2.8463]
%           std_bHPE: [0.1008 0.3329 0.9010]
% 
% resadd = 3
% 
% resadd =
% 
%      3
% 
% 
% ans = 
% 
%   struct with fields:
% 
%             theory: [2.4977 17.6753 7.0693]
%              convE: [0.7385 12.1046 1.4692]
%          std_convE: [0.0148 0.1873 0.0310]
%         L1_locstat: NaN
%     std_L1_locstat: NaN
%                HPE: [0.7853 11.8981 1.7900]
%            std_HPE: [0.0376 0.1851 0.4965]
%               bHPE: [0.7841 12.1355 2.2171]
%           std_bHPE: [0.1135 0.3288 1.0025]
% 
% 
% ans = 
% 
%   struct with fields:
% 
%             theory: [6.2442 12.9965 9.0169]
%              convE: [5.5809 10.5171 6.4517]
%          std_convE: [0.0757 0.1886 0.2483]
%         L1_locstat: 5.5776
%     std_L1_locstat: 0.4942
%                HPE: [5.5980 10.2809 5.6518]
%            std_HPE: [0.0597 0.2190 0.4811]
%               bHPE: [5.6222 10.4451 6.3282]
%           std_bHPE: [0.0893 0.4189 0.8415]