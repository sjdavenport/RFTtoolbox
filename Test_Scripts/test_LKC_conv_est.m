%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%    This script tests the LKC_conv_est.m function
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all

addpath( genpath( "/home/drtea/matlabToolboxes/RFTtoolbox/" ) )


%% %-----------------------------------------------------------------------
%    Test section D = 1
%--------------------------------------------------------------------------
% parameter for the field
T      = 50;
nsubj  = 500;
FWHM   = sigma2FWHM(20);
resAdd = 7;
pad    = ceil( 4*FWHM2sigma( FWHM ) );
% method = "numerical";
method = "analytical";

mask = zeros( [ T+2*pad 1 ] );
mask( (pad+1):(end-pad) ) = 1;

% get true LKC from continuous theory
trueL = LKC_isogauss_theory( FWHM, T  );

% generate test data
lat_data = randn( [ T+2*pad nsubj ] );

% closest approximation of the continuous field
LKC = LKC_conv_est( lat_data, mask, FWHM, resAdd, 0,...
                    ceil( resAdd / 2 ), method );
[ trueL; LKC.hatL ]


% estimate of LKCs
LKC1 = LKC_conv_est( lat_data, mask, FWHM, resAdd, 0, 0 )
LKC2 = LKC_conv_est( lat_data, mask, FWHM, resAdd, 0, 1 )
LKC3 = LKC_conv_est( lat_data, mask, FWHM, resAdd, 1 )


%% %-----------------------------------------------------------------------
%    Test section D = 2
%--------------------------------------------------------------------------
% parameter for the field
T      = 49;
nsubj  = 100;
FWHM   = sigma2FWHM(5);
resAdd = 1;
mask_opt = [1 1];
pad = ceil( 4*FWHM2sigma( FWHM ) );

mask = zeros( [ T+2*pad T+2*pad ] );
mask( (pad+1):(end-pad), (pad+1):(end-pad) ) = 1;

% get true LKC
trueL = LKC_isogauss_theory( FWHM, [ T T ] );

% generate test data
lat_data = randn( [ T+2*pad T+2*pad nsubj ] );

% closest approximation of the continuous field
LKC = LKC_conv_est( lat_data, mask, FWHM, resAdd, 0 );
[ trueL; LKC.hatL ]

% estimate of LKCs using different options for the mask
LKC1 = LKC_conv_est( lat_data, mask, FWHM, resAdd, 1 )
LKC2 = LKC_conv_est( lat_data, mask, FWHM, resAdd, 0, 0 )
LKC3 = LKC_conv_est( lat_data, mask, FWHM, resAdd, 0, 1 )


%% %-----------------------------------------------------------------------
%    Test section D = 3
%--------------------------------------------------------------------------
% parameter for the field
T      = 30;
nsubj  = 50;
FWHM   = sigma2FWHM(5);
resAdd = 1;
pad    = ceil( 4*FWHM2sigma( FWHM ) );

% generate rectangular mask with a padded zero collar 
mask = zeros( [ T+2*pad T+2*pad T+2*pad ] );
mask( (pad+1):(end-pad), (pad+1):(end-pad), (pad+1):(end-pad) ) = 1;

% get true LKC
trueL = LKC_isogauss_theory( FWHM, [ T T T ] );

% generate test data
lat_data = randn( [ T+2*pad T+2*pad T+2*pad nsubj ] );

% closest approximation of the continuous field
tic
LKC = LKC_conv_est( lat_data, mask, FWHM, resAdd, 0 );
toc
[ trueL; LKC.hatL ]

% estimate of LKCs
LKC1 = LKC_conv_est( lat_data, mask, FWHM, resAdd, 1, 3 )
LKC2 = LKC_conv_est( lat_data, mask, FWHM, resAdd, 0, 0 )
LKC3 = LKC_conv_est( lat_data, mask, FWHM, resAdd, 1 )
