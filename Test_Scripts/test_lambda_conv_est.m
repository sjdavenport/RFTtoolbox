%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%    This script tests the Lambda_conv_est.m function
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all

addpath( genpath( "/home/drtea/matlabToolboxes/RFTtoolbox/" ) )

%% Test section D=1
%--------------------------------------------------------------------------

%%% show that the Gaussian kernel does the same as using the naive method
% without asuming seperable kernels
% Field parameters
FWHM  = 10;
D     = 1;
Nsubj = 120;
T    = 100;
Dim  = ones( [ 1 D ] ) * T;
siz  = ceil( 4* FWHM2sigma( FWHM ) );
mask = ones( [ Dim 1 ] );
sM   = size( mask );

% resolution parameters
resAdd = 1;

% generate data
Y = randn( [ Dim(1) Nsubj ] );

%%%% Lambda est function
enlarge = 0;
h = 1e-5;
hatLambda_conv = Lambda_conv_est( Y, FWHM, resAdd, enlarge );
%hatLambda_num  = Lambda_numeric_est( Y, FWHM, resAdd, enlarge, h );
%Lambda_array   = Lambda_est( Y, FWHM, D, resAdd, h );

%sameArray( hatLambda_conv, hatLambda_num )
%sameArray( hatLambda_conv, Lambda_array )
%sameArray( Lambda_array, hatLambda_num )

%% Test section D=2
%--------------------------------------------------------------------------

%%%%%% show that the Gaussian kernel does the same as using the naive method
%%% without asuming seperable kernels
% Field parameters
FWHM  = 10;
D     = 2;
Nsubj = 120;
T     = 100;
Dim   = ones( [ 1 D ] ) * T;
siz   = ceil( 4* FWHM2sigma( FWHM ) );
mask  = ones( Dim );
sM    = size(mask);

% resolution parameters
resAdd = 1;
dx = 1 / ( resAdd + 1 );
Dimhr = ( sM - 1 ) * resAdd + sM;

% generate data
Y = randn( [ Dim Nsubj ]);

%%%% Lambda est function
Lambda_est = Lambda_conv_est( Y, FWHM, resAdd, 3 );

size(Lambda_est)
Dimhr


%% Test section D=3
%--------------------------------------------------------------------------

%%%%%% show that the Gaussian kernel does the same as using the naive
%%%%%% method without asuming seperable kernels
% Field parameters
FWHM  = 3;
D     = 3;
Nsubj = 120;
T     = 40;
Dim   = ones( [ 1 D ] ) * T;
siz   = ceil( 4* FWHM2sigma( FWHM ) );
mask  = ones( Dim );
sM    = size( mask );

% resolution parameters
resAdd = 1;
dx = 1 / ( resAdd + 1 );
Dimhr = ( sM - 1 ) * resAdd + sM;

% generate data
Y = randn( [ Dim Nsubj ]);

%%%% Lambda est function
Lambda_est = Lambda_conv_est( Y, FWHM, resAdd, 6 );

size(Lambda_est)
Dimhr