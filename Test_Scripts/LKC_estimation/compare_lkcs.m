%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script compares LKC estimation using different methods
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1D
FWHM = 6; resadd = 1; nsubj = 50; nvox = 100;
lat_data = wnfield(nvox, nsubj);

% Convolution L_1 (voxmndest)
cfield  = convfield_Field( lat_data, FWHM, 0, resadd );
dcfield = convfield_Field( lat_data, FWHM, 1, resadd );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );

% HPE L_1
HPE  = LKC_HP_est( cfield, 1, 1 );

% Theory
L_theory = LKC_isogauss_theory( FWHM, nvox );

% SPM
L_spm = LKC_SPM_est( FWHM, lat_data.mask );

% Compare results
struct( 'theory', L_theory, 'voxmfd_est', L_conv, ...
        'HPE', HPE.hatL, 'SPM', L_spm )
    
% FT: Works perfectly well, now.

%% 2D (no padding) big mismatch between the estimators here
% Note: can reduce resadd for fast implementation of course)
FWHM = 6; nsubj = 20; Dim = [49,49];
resadd = 9; lat_data = wnfield(Dim, nsubj);

% Convolution L_1 (voxmndest)
cfield  = convfield_Field( lat_data, FWHM, 0, resadd );
dcfield = convfield_Field( lat_data, FWHM, 1, resadd );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );


% Convolution L_1 (LKC_conv_est), not the same as above for some reason
% FT: this function because I did not developed further had a bug, which
% was a bug in mask_highres, which gave all 1's weights even if resolution
% is increased, which is simply wrong.
LKCs = LKC_conv_est( lat_data.field, lat_data.mask, FWHM, resadd );
LKCs.hatL;

% Convolution (previous implemenation) gives a very different answer
Gker_param = FWHM2sigma(FWHM); D = length(Dim);
L_oldconv = LKC_GaussConv( lat_data.field, Gker_param, D, resadd );

% HPE seems to match the LKC_GaussConv answer here
HPE  = LKC_HP_est( cfield, 1, 1);
newHPE = HPE.hatL';

% In this case old HPE matches new HPE
D = length(Dim);
HPE  = LKCestim_HPE( cfield.field, D, cfield.mask, 1);
oldHPE = HPE.hatn';

% SPM (off of course because of non-stationarity but certainly more similar
% to the HPE and LKCGaussconv estimators)
[ L_spm, L0 ] = LKC_SPM_est( FWHM, lat_data.mask );

% IsoGauss Theory
L_theory = LKC_isogauss_theory( FWHM, Dim );

% Compare results
struct( 'theory', L_theory, 'voxmfd_est', L_conv, 'conv_est', LKCs.hatL,...
        'old_conv', L_oldconv ,...
        'newHP', newHPE,...
        'oldHPE', oldHPE, 'spm', L_spm )

    % FT: as seen later HP might have a small downward bias, that's why
    
    
%% 2D (with padding to ensure stationarity)
FWHM = 12; nsubj = 100; Dim = [49,49]; D = length(Dim);
resadd = 5; pad = ceil( 4 * FWHM2sigma( FWHM ) ); mask = ones(Dim);
padded_mask = logical( pad_vals( mask, pad) ); lat_data = wnfield(padded_mask, nsubj);

% Convolution L_1 (voxmndest) seems to work in this case
lat_masked = 0;
cfield  = convfield_Field( lat_data, FWHM, 0, resadd, lat_masked );
dcfield = convfield_Field( lat_data, FWHM, 1, resadd, lat_masked );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );

%Note LKC_conv_est can't be used in this case as there is no way to
%ensure the lat_masked = 0 option.

% HPE (This seems to be off, I'm not sure what is wrong here?)
HPE  = LKC_HP_est( cfield, 1, 1);
newHPEestimate = HPE.hatL';

% HPE old version
HPE  = LKCestim_HPE( cfield.field, D, cfield.mask, 1);
oldHPEestimate = HPE.hatn';

% SPM
[ L_spm, L0 ] = LKC_SPM_est( FWHM, mask );

% IsoGauss Theory
L_theory = LKC_isogauss_theory( FWHM, Dim );

% Compare results
struct( 'theory', L_theory, 'voxmfd_est', L_conv, 'newHP', newHPEestimate,...
        'oldHPE', oldHPEestimate, 'spm', L_spm )
    
% FT again seems to work well for convoltuiton estimator and large FHWM things seem to converge.

%% 3D (no padding) box
FWHM = 10; nsubj = 20; Dim = [5,5,5];
resadd = 1; mask = true(Dim); lat_data = wnfield(mask, nsubj);

% Convolution L_1 (voxmndest)
cfield  = convfield_Field( lat_data, FWHM, 0, resadd );
dcfield = convfield_Field( lat_data, FWHM, 1, resadd );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );

% HPE (Worryingly L_3 can be negative here sometimes)
HPE  = LKC_HP_est( cfield, 1, 1);
newHPEestimate = HPE.hatL';

% isogauss theory (Off of course due to the edge correction)
L_theory = LKC_isogauss_theory( FWHM, Dim );

% SPM (Off of course due to the edge correction)
[ L_spm, L0 ] = LKC_SPM_est( FWHM, mask );

% Compare results
struct( 'theory', L_theory, 'voxmfd_est', L_conv, 'newHP', newHPEestimate,...
        'spm', L_spm )

%% 3D (with padding) box
FWHM = 3; nsubj = 20; Dim = [5,5,5]; D = length(Dim);
resadd = 3; pad = ceil( 4 * FWHM2sigma( FWHM ) ); mask = ones(Dim);
padded_mask = logical( pad_vals( mask, pad) ); lat_data = wnfield(padded_mask, nsubj);

% Convolution L_1 (voxmndest) seems to work in this case
lat_masked = 0;
enlarge = ceil(resadd/2);
% enlarge = 0;
cfield  = convfield_Field( lat_data, FWHM, 0, resadd, lat_masked, enlarge );
dcfield = convfield_Field( lat_data, FWHM, 1, resadd, lat_masked, enlarge );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );

% HPE (This seems to be off, I'm not sure what is wrong here?)
HPEnew  = LKC_HP_est( cfield, 1, 1);

% HPE old version
HPEold  = LKCestim_HPE( cfield.field, D, cfield.mask, 1);

% isogauss theory
L_theory = LKC_isogauss_theory( FWHM, Dim );

% SPM theory (need to set enlarge = 0 for SPM to be correct here)
[ L_spm, L0 ] = LKC_SPM_est( FWHM, mask );

% Compare results
struct( 'theory', L_theory, 'voxmfd_est', L_conv, 'newHP', HPEnew.hatL' ,...
        'oldHPE', HPEold.hatn', 'spm', L_spm )

%% Compare theory estimates, different because SPM is coded for enlarge = 0 it seems
Dim = [5,5,5];

% isogauss theory
LKC_isogauss_theory( FWHM, Dim )

% isogauss theory enlarge = 0, or equivalently add +1 in the SPM theory
LKC_isogauss_theory( FWHM, Dim-1 )

% SPM theory (incorrect)
[ L_spm, L0 ] = LKC_SPM_est( FWHM, ones(Dim) );
L_spm

%% SPM almost matches for large domains (as we would expect)
Dim = [100,100,100];

% isogauss theory
LKC_isogauss_theory( FWHM, Dim )

% SPM theory
[ L_spm, L0 ] = LKC_SPM_est( FWHM, ones(Dim) );
L_spm

%% 3D boundary sphere example (just Conv, HPE)
% FT: I fixed L1, but I think it is still not completely correct
T = 5;
D = 3;
Dim = T * ones( [ 1 D ] );
FWHM = 5;
resadd = 3;
nsubj = 200; % FT:for pad= 0 change nsubj to see that negative values appear more often for small sample sizes
pad   = 0; %ceil( 4 * FWHM2sigma( FWHM ) ); %  for non-stationarity 

% Change to false to see that everything works well in stationary approach
lat_mask = false;

% Cube manifold
mask = true( Dim );
mask = logical( pad_vals( mask, pad) );
cubeL = LKC_isogauss_theory( FWHM, Dim );

% Sphere manifold
mask_sphere = true( Dim );
mask_sphere = bndry_voxels( logical( pad_vals( mask_sphere, pad) ), 'full' );
sm = size( mask_sphere );

lambda   = 1 / ( 2 * FWHM2sigma( FWHM )^2 );
volmask  = 2*T^2 + 2*(T-2)*T + 2* (T-2)^2;
surfmask = 6 * ( T^2 + (T-2)^2 );
linemask = 12 * T / 4 - 12 * ( T - 2 ) * 1 / 4;

sphere = VoxManifold( mask_highres( mask_sphere, 1 ), ...
                      xvals_highres( {1:sm(1),1:sm(2),1:sm(3)}, 1 ) );
sphereL  = [ sqrt(lambda), lambda, lambda^(3/2) ] .* LKC_est( sphere );
sphereL1 = [ sqrt(lambda), lambda, lambda^(3/2) ] .* [ NaN surfmask/2 volmask ];

% Get data
lat_data = wnfield( mask_sphere, nsubj );

% Convolution L_1 (voxmndest)
cfield  = Mask( convfield_Field( lat_data, FWHM, 0, resadd, lat_mask ) );
dcfield = Mask( convfield_Field( lat_data, FWHM, 1, resadd, lat_mask ) );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );

% Estimate using christoffel symbols
d2cfield = Mask( convfield_Field( lat_data, FWHM, 2, resadd, false ) );
[L_conv_2,L0_conv] = LKC_voxmfd_est( cfield, dcfield, d2cfield );

% HPE old (and correct) version
HPEold  = LKCestim_HPE( cfield.field, D, cfield.mask, 1 );

% HPE (I'm assuming this is off based on previous examples, but have included for completeness)
HPEnew  = LKC_HP_est( cfield, 1, 1);

% Compare results
struct( 'theory', sphereL, 'voxmfd_est', L_conv, 'newHP', HPEnew.hatL' ,...
        'oldHPE', HPEold.hatn' )

%% Calculate average estimates
niters = 50;
Dim = [5,5,5]; FWHM = 3; resadd = 9;
mask_sphere = true( Dim ); D = length(Dim);
mask_sphere = bndry_voxels( mask_sphere, 'full');

L_conv_average = zeros(1,D);
L_conv_2_average = zeros(1,D);
L_HPE_average = zeros(1,D);

% Takes a few minutes (to reduce the time can decrease resadd)
for I = 1:niters
    I
    lat_data = wnfield(mask_sphere, nsubj);
    
    % Convolution L_1 (voxmndest)
    cfield  = convfield_Field( lat_data, FWHM, 0, resadd );
    dcfield = convfield_Field( lat_data, FWHM, 1, resadd );
    [L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield );
    L_conv_average = L_conv_average + L_conv;
    
    dcfield2 = convfield_Field( lat_data, FWHM, 2, resadd );
    [L_conv_2,L0_conv] = LKC_voxmfd_est( cfield, dcfield, dcfield2 );
    L_conv_2_average = L_conv_2_average + L_conv_2;
    
    % HPE old (and correct) version
    HPE  = LKCestim_HPE( cfield.field, D, cfield.mask, 1);
    L_HPE_average = HPE.hatn' + L_HPE_average;
end

L_conv_average/niters
L_conv_2_average/niters
L_HPE_average/niters

% When I run the average loop I got
% L_conv_average/niters = 21.3467   23.8655    6.7526  and 
% L_conv_2_average/niters =  21.4748   23.8655    6.7526  and 
% L_HPE_average/niters = 2.0199   22.2506    5.5273

% This indicates that L_2 and L_3 are likely computed correctly while L_1 is off
% in the convolution estimator.

%% 3D MNImask LKCs
FWHM = 3; nsubj = 20; resadd = 1; 
mask = getMNImask; lat_data = wnfield(mask, nsubj);

% Convolution L_1 (voxmndest)
cfield  = convfield_Field( lat_data, FWHM, 0, resadd );
dcfield = convfield_Field( lat_data, FWHM, 1, resadd );
[L_conv,L0_conv] = LKC_voxmfd_est( cfield, dcfield )

% HPE (Worryingly L_3 can be negative here sometimes)
HPE  = LKC_HP_est( cfield, 1, 1);
newHPEestimate = HPE.hatL'

% SPM (Off of course due to the edge correction)
[ L_spm, L0 ] = LKC_SPM_est( FWHM, mask )

%% 3D MNImask LKCs with padding to ensure stationarity
FWHM = 3; nsubj = 20; mask = getMNImask;
pad = ceil( 4 * FWHM2sigma( FWHM ) ); padded_mask = logical( pad_vals( mask, pad) );
resadd = 1; lat_data = wnfield(padded_mask, nsubj);
lat_masked = 0;

% Convolution L_1 (voxmndest)
[tcfield, cfields] = convfield_t_Field( lat_data, FWHM, resadd, lat_masked );
dcfield = convfield_Field( lat_data, FWHM, 1, resadd, lat_masked );
[L_conv,L0_conv] = LKC_voxmfd_est( cfields, dcfield )

% HPE (Worryingly L_3 can be negative here sometimes)
HPE  = LKC_HP_est( cfields, 1, 1);
newHPEestimate = HPE.hatL'
 
% SPM (Off of course due to the edge correction)
[ L_spm, L0 ] = LKC_SPM_est( FWHM, mask )

%% Compare EC curves
[ curve, thresholds ] = ECcurve( tcfield, [-5,5]);

EEC_spm = EEC_calc( thresholds, L_spm, L0, 'T', nsubj )
EEC_hpe = EEC_calc( thresholds, newHPEestimate, L0, 'T', nsubj )
EEC_conv = EEC_calc( thresholds, L_conv, L0, 'T', nsubj )
EEC_mix = EEC_calc( thresholds, [newHPEestimate(1), L_conv(2), L_conv(3)], L0, 'T', nsubj )

plot(thresholds, curve)
hold on
plot(thresholds, EEC_spm)
hold on 
plot(thresholds, EEC_hpe)
plot(thresholds, EEC_conv)
plot(thresholds, EEC_mix)
legend('Observed', 'SPM', 'HPE', 'Conv', 'Mix')

%%
Lambda = Lambda_theory( FWHM, 3 );
sum(MNImask(:))*prod(diag(Lambda))^(1/2)

%% MNI SPM highres
mask = logical(imgload('MNImask')); resadd = 3;
MNIhr = mask_highres(mask, resadd);

[ L_spm, L0 ] = LKC_SPM_est( FWHM, MNIhr)

%% MNI SPM highres
mask = logical(imgload('MNImask')); resadd = 1;
MNIhr = mask_highres(mask, resadd);

[ L_spm, L0 ] = LKC_SPM_est( FWHM, MNIhr )
