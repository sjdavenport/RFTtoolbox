%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the Lambda_numeric_est function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% 1D Examples
%% 1D comparison with stationarity
FWHM = 10; D = 1; sigma = FWHM2sigma(FWHM); nvox = 100; nsubj = 10000;
Lambda_theory = diag(repmat(sigma^(-2),1,D))/2; lat_data = normrnd(0,1,nvox,nsubj);
Lambda_estimates = Lambda_numeric_est( lat_data, FWHM, D );
Lambda_theory
mean(Lambda_estimates(5:95)) % Not the ends because of edge effects
plot(Lambda_estimates); hold on; abline('h', Lambda_theory);

%% 2D (stationary example)
FWHM = 5; sigma = FWHM2sigma(FWHM); Dim = [100,100]; D = length(Dim); nsubj = 50;
Lambda_theory = diag(repmat(sigma^(-2),1,D))/2; lat_data = normrnd(0,1,[Dim,nsubj]);
Lambda_estimates = Lambda_numeric_est( lat_data, FWHM, D );
first_entry = Lambda_estimates(5:95,5:95,1,1);
Lambda_theory(1,1)
mean(first_entry(:))
onetwo_entry = Lambda_estimates(5:95,5:95,1,2);
Lambda_theory(1,2)
mean(onetwo_entry(:))
plot(first_entry(:)); hold on; abline('h', Lambda_theory(1,1));

%% 3D (stationary example) (takes around a minute)
FWHM = 3; sigma = FWHM2sigma(FWHM); Dim = [50,50,50]; D = length(Dim); nsubj = 50;
Lambda_theory = diag(repmat(sigma^(-2),1,D))/2; lat_data = normrnd(0,1,[Dim,nsubj]);
resadd = 0;
%%
tic; Lambdahat_num = Lambda_numeric_est( lat_data, FWHM, 0 ); toc
tic; Lambdahat_conv = Lambda_conv_est( lat_data, FWHM, 0 ); toc
[is_same, array_diff] = same_array(Lambdahat_conv, Lambdahat_num)

%%
first_entry = Lambdahat_num(5:45,5:45,5:45,1,1);
Lambda_theory(1,1)
mean(first_entry(:))
onetwo_entry = Lambdahat_num(5:45,5:45,5:45,1,2);
Lambda_theory(1,2)
mean(onetwo_entry(:))
plot(first_entry(:)); hold on; abline('h', Lambda_theory(1,1));

%%
[is_same, array_diff] = same_array(Lambdahat_conv, Lambdahat_num)
first_entry = Lambdahat_conv(5:45,5:45,5:45,1,1);
Lambda_theory(1,1)
mean(first_entry(:))
onetwo_entry = Lambdahat_conv(5:45,5:45,5:45,1,2);
Lambda_theory(1,2)
mean(onetwo_entry(:))

%% Visual comparison to Lambda_conv_est
FWHM = 3; sample_size = 50; Dim = [5,5]; mask = true(Dim); resadd = 11;
nsubj = 1000; lat_data = normrnd(0,1,[Dim, nsubj]);
a = 2; b = 2;
ghat = Lambda_numeric_est( lat_data, FWHM, resadd );
first_entry = ghat(:,:,a,b);
mean(ghat(:))
subplot(2,1,1);imagesc(first_entry); title('Lambda num est')
ghat = Lambda_conv_est( lat_data, FWHM, resadd );
first_entry = ghat(:,:,a,b);
mean(ghat(:))
subplot(2,1,2); imagesc(first_entry); title('Lambda conv est')

%% %% Comparison to Lambda_conv_est
%% 1D example
FWHM  = 10;
D     = 1;
Nsubj = 120;
T    = 100;
Dim  = ones( [ 1 D ] ) * T;
siz  = ceil( 4* FWHM2sigma( FWHM ) );
mask = ones( [ Dim 1 ] );
sM   = size( mask );

% resolution parameters
resadd = 1;

% generate data
Y = randn( [ Dim(1) Nsubj ] );

%%%% Lambda est function
enlarge = 0;

tic
Lambdahat_conv = Lambda_conv_est( Y, FWHM, resadd, enlarge );
toc
tic
Lambdahat_num = Lambda_numeric_est( Y, FWHM, resadd, enlarge );
toc
plot(Lambdahat_num, Lambdahat_conv)
same_array(Lambdahat_num, Lambdahat_conv)

%% 2D example
FWHM  = 10;
D     = 2;
Nsubj = 120;
T     = 100;
Dim   = ones( [ 1 D ] ) * T;
siz   = ceil( 4* FWHM2sigma( FWHM ) );
mask  = ones( Dim );
sM    = size(mask);

% resolution parameters
resdd = 1;
dx = 1 / ( resadd + 1 );
Dimhr = ( sM - 1 ) * resadd + sM;

% generate data
Y = randn( [ Dim Nsubj ]);

%%%% Lambda est function
Lambdahat_conv = Lambda_conv_est( Y, FWHM, resadd, 3 );
Lambdahat_num = Lambda_numeric_est( Y, FWHM, resadd, 3 );

size(Lambdahat_conv)
size(Lambdahat_num)
[is_same, array_diff] = same_array(Lambdahat_conv, Lambdahat_num)

%% 3D example
% Field parameters
FWHM  = 3;
D     = 3;
Nsubj = 40;
T     = 40;
Dim   = ones( [ 1 D ] ) * T;
siz   = ceil( 4* FWHM2sigma( FWHM ) );
mask  = ones( Dim );
sM    = size( mask );

% resolution parameters
resadd = 0;
dx = 1 / ( resadd + 1 );
Dimhr = ( sM - 1 ) * resadd + sM;

% generate data
Y = randn( [ Dim Nsubj ]);

%%%% Lambda est function
tic
Lambdahat_conv = Lambda_conv_est( Y, FWHM, resadd, 6 );
toc
tic
Lambdahat_num = Lambda_numeric_est( Y, FWHM, resadd, 6 );
toc

size(Lambdahat_conv)
size(Lambdahat_num)
[is_same, array_diff] = same_array(Lambdahat_conv, Lambdahat_num)

