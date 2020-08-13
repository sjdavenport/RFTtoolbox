%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the record_coverage when the fields are not
%%%    Gaussian
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D examples
%% 1D white t example
FWHM = 3; sample_size = 100; nvox = 100; resadd = 3; df = 3;
spfn = @(nsubj) wtfield( nvox, nsubj, df ); niters = 1000;
params = ConvFieldParams( FWHM, resadd );
coverage = record_coverage( spfn, sample_size, params, niters)
% ECcurveanal(coverage, true(nvox,1), sample_size, 0.1)

%% 1D bootstrap t example
rng(1)
nvox = 100; D = 1;
spfn = @(nsubj) wtfield(nvox, nsubj, 3);
data = spfn(560).field;
spfn = get_sample_fields( data, true(1,nvox)', D );

FWHM = 3; resadd = 1; sample_size = 200;  niters = 1000;
params = ConvFieldParams( FWHM, resadd );
coverage = record_coverage( spfn, sample_size, params, niters)

%% 1D white Laplacian example
FWHM = 3; sample_size = 20; nvox = 20; resadd = 1; scale = 2;
spfn = @(nsubj) wlfield( nvox, nsubj, scale ); niters = 1000;
params = ConvFieldParams( FWHM, resadd );
coverage = record_coverage( spfn, sample_size, params, niters)
% ECcurveanal(coverage, true(nvox,1), sample_size, 0.1)

%% Second small 1D example
FWHM = 3; sample_size = 10; nvox = 10; resadd = 1;
spfn = @(nsubj) wnfield( nvox, nsubj ); niters = 1000;
params = ConvFieldParams( FWHM, resadd );
record_coverage( spfn, sample_size, params, niters)

%% 1D bootstrap example
nvox = 100; tot_nsubj = 600; D = 1;
data = normrnd(0,1,nvox, tot_nsubj);
spfn = get_sample_fields( data, true(1,nvox)', D );
FWHM = 3; sample_size = 20; resadd = 1; niters = 1000;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
record_coverage( spfn, sample_size, params, niters)

%% %% 2D examples
%% small 2D white t example
%(Note that on a small domain conv and finelat are similar even for reasonable resadd)
Dim = [10,10]; FWHM = 3; resadd = 1; sample_size = 100;
spfn = @(nsubj) wtfield(Dim, nsubj, 3); niters = 1000;
params = ConvFieldParams( [FWHM, FWHM], resadd );
coverage = record_coverage( spfn, sample_size, params, niters)
% ECcurveanal(coverage, true(Dim), sample_size, 0.1)

%% 2D t bootstrap example
%% Set Sample field function and parameters
D = 2; niters = 1000;
spfn = get_sample_fields( im_store_2D, MNImask_2D, D );
FWHM = 3; resadd = 1;

Dim = [10,10]; 
spfn = @(nsubj) wtfield(Dim, nsubj, 3); niters = 1000;
data = spfn(560).lat_data.field;

FWHM = 3; resadd = 1; sample_size = 200;
params = ConvFieldParams( [FWHM, FWHM], resadd );
coverage = record_coverage( spfn, sample_size, params, niters)
% ECcurveanal(coverage, true(Dim), sample_size, 0.1)

%% small 2D white Laplacian example
%(Note that on a small domain conv and finelat are similar even for reasonable resadd)
Dim = [30,30]; FWHM = 3; resadd = 1; scale = 0.01; sample_size = 100;
spfn = @(nsubj) wlfield(Dim, nsubj, scale); niters = 1000;
params = ConvFieldParams( [FWHM, FWHM], resadd );
coverage = record_coverage( spfn, sample_size, params, niters)
% ECcurveanal(coverage, true(Dim), sample_size, 0.1)

%% small 2D example - using HPE
FWHM = 3; resadd =  11; sample_size = 50;
spfn = @(nsubj) wnfield(Dim, nsubj); niters = 1000;
record_coverage_dep( spfn, sample_size, FWHM, resadd, niters, 'hpe')

%% 2D example on MNImask
MNImask = imgload('MNImask');
MNImask_2D = MNImask(:,:,45);
FWHM = 3; resadd =  1; sample_size = 20;
params = ConvFieldParams( [FWHM,FWHM], resadd );
spfn = @(nsubj) wnfield(logical(MNImask_2D), nsubj); niters = 1000;
coverage = record_coverage( spfn, sample_size, params, niters)
ECcurveanal(coverage, MNImask_2D, sample_size, 0.5)

%% 2D bootstrap example on MNImask
MNImask = imgload('MNImask');
MNImask_2D = MNImask(:,:,45);
FWHM = 3; resadd =  1; sample_size = 20;
params = ConvFieldParams( [FWHM,FWHM], resadd );
spfn = @(nsubj) wnfield(logical(MNImask_2D), nsubj); niters = 1000;
coverage = record_coverage( spfn, sample_size, params, niters)
ECcurveanal(coverage, MNImask_2D, sample_size, 0.5)


%% 2D bootstrap
Dim = [30,30]; tot_nsubj = 600; D = length(Dim);
data = normrnd(0,1,[Dim, tot_nsubj]);
spfn = get_sample_fields( data, true(Dim), D );
FWHM = 3; sample_size = 20; resadd = 1; niters = 1000;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
record_coverage( spfn, sample_size, params, niters)

%% 2D bootstrap (high sample to tot subj ratio)
Dim = [40,40]; tot_nsubj = 100; D = length(Dim);
data = normrnd(0,1,[Dim, tot_nsubj]);
spfn = get_sample_fields( data, true(Dim), D );
FWHM = 3; sample_size = 20; resadd = 1; niters = 1000;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
record_coverage( spfn, sample_size, params, niters)

%% large 2D example 
%(Note that on a small domain conv and finelat are similar even for reasonable resadd)
FWHM = 3; resadd = 1; sample_size = 50; Dim = [50,50];
spfn = @(nsubj) wnfield(Dim, nsubj); niters = 1000;
params = ConvFieldParams( [FWHM,FWHM], resadd );
record_coverage( spfn, sample_size, params, niters)

%% small 2D sphere example
FWHM = 3; resadd = 3; sample_size = 50; Dim = [5,5];
mask = bndry_voxels(true(Dim), 'full');
spfn = @(nsubj) wnfield(mask, nsubj); niters = 1000;
params = ConvFieldParams( [FWHM,FWHM], resadd );
coverage = record_coverage( spfn, sample_size, params, niters);

%% Compare max distribution and EC curves
L0 = EulerChar(mask, 0.5, 2); %@Fabian this is wrong
L0 = 1;
thresholds = -6:0.05:6;
av_L_est = mean(coverage.storeLKCs,2)';
EEC_conv = EEC( thresholds, av_L_est, L0, 'T', sample_size -1 );
plot(EEC_conv)
%% large 2D example
FWHM = 3; sample_size = 50; Dim = [50,50]; resadd = 1;
spfn = @(nsubj) wnfield(Dim, nsubj); niters = 1000;
record_coverage( spfn, sample_size, FWHM, resadd, niters)
 
%% %% 3D examples
%% Small 3D example
FWHM = 3; sample_size = 50; Dim = [5,5,5]; resadd = 1;
spfn = @(nsubj) wnfield(Dim, nsubj); niters = 1000;
record_coverage( spfn, sample_size, FWHM, resadd, niters)

%% Large 3D example (takes a while)
FWHM = 3; sample_size = 50; Dim = [50,50,50]; resadd = 1;
spfn = @(nsubj) wnfield(Dim, nsubj); niters = 1000;
record_coverage( spfn, sample_size, FWHM, resadd, niters, 'conv', 0)

%% MNImask example
mask = logical(imgload('MNImask')); FWHM = 3; sample_size = 10; 
Dim = [91,109,91]; resadd = 1; spfn = @(nsubj) wnfield(Dim, nsubj); 
niters = 1000;
record_coverage( spfn, sample_size, FWHM, resadd, niters, 'conv', 0)

%% Non-stat 1D
dim = [50 1]; D = 1;
voxmap = randsample(50,50,0);
spfn = @(nsubj) cnfield( dim, 10, voxmap, 0, nsubj );

FWHM = 3; sample_size = 20; resadd = 3; niters = 1000;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
cov_conv = record_coverage( spfn, sample_size, params, niters)

%% Non-stat 1D HPE
dim = [50 1]; D = 1;
voxmap = randsample(50,50,0);
spfn = @(nsubj) cnfield( dim, 10, voxmap, 0, nsubj );

FWHM = 3; sample_size = 20; resadd = 3; niters = 1000;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
cov_HPE = record_coverage( spfn, sample_size, params, niters, 'HPE')

%% Non-stat 1D
dim = [50 1]; D = 1;
voxmap = randsample(50,50,0);
spfn = @(nsubj) cnfield( dim, 10, voxmap, 0, nsubj );

FWHM = 3; sample_size = 20; resadd = 3; niters = 1000;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
coverage1D = record_coverage( spfn, sample_size, params, niters)

%% Examine EC curves
LKC_estimate = mean(coverage1D.storeLKCs,2)';
L0 = 1;
[ curve, x ] = maxECcurve( coverage1D.convmaxima, 0.1 )
[ curve_all, x_all ] = maxECcurve( coverage1D.allmaxima, 0.1 )

plot(x, curve);
curve_conv = EEC( x, LKC_estimate, L0, 'T', sample_size -1 );
hold on 
plot(x,curve_conv)
hold on 
plot(x_all,curve_all)
legend('max', 'LKC', 'all')


%% Examine EC curves
LKC_estimate = mean(coverage1D.storeLKCs,2)';
L0 = 1;
[ curve, x ] = maxECcurve( coverage1D.allmaxima, 0.1 )
plot(x, curve);
curve_conv = EEC( x, LKC_estimate, L0, 'T', sample_size -1 );
hold on 
plot(x,curve_conv)

%%


%% Examine tstat
dim = [100 1]; D = 1;
voxmap = 1:dim;
tmp2 = mod( voxmap, 6);
voxmap = [ voxmap( tmp2 < 2 ), voxmap( tmp2 >= 2 ) ];
spfn = @(nsubj) cnfield( dim, 10, voxmap, 0, nsubj );

params = ConvFieldParams( repmat(FWHM,1,D), resadd );
tcfield = convfield_t( spfn(100), params );
plot(tcfield)

%% Non-stat 2D
Dim = [50 50]; D = length(Dim);
voxmap = randsample(prod(Dim),prod(Dim),0);
spfn = @(nsubj) cnfield( Dim, 10, voxmap, 0, nsubj );

FWHM = 3; sample_size = 100; resadd = 1; niters = 1000;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
record_coverage( spfn, sample_size, params, niters)

%% Gaussian multiplier fields
Dim = [50,1]; nsubj = 20; D = 1;
initial_data = wnfield(Dim, nsubj);
spfn = @(nsubj) Gmult(initial_data); %Gmult generates a Gaussian multiplier field
FWHM = 3; sample_size = 20; resadd = 1; niters = 1000;
params = ConvFieldParams( repmat(FWHM,1,D), resadd );
coverage1D = record_coverage( spfn, sample_size, params, niters)

%% EC curve analysis
LKC_estimate = mean(coverage1D.storeLKCs,2)';
L0 = EulerChar(mask_slice, 0.5, D);
[ curve, x ] = maxECcurve( coverage1D.convmaxima, 0.1 )
[ curve_all, x_all ] = maxECcurve( coverage1D.allmaxima, 0.1 )

plot(x, curve);
curve_conv = EEC( x, LKC_estimate, L0, 'T', sample_size -1 );
hold on 
plot(x,curve_conv)
hold on 
plot(x_all,curve_all)
legend('max', 'LKC', 'all')