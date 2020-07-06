clear all
close all

% General parameters
sigma = 5;
T     = 49;
FWHM  = sigma2FWHM( 5 );
nsubj = 20;
pad   = ceil(4*sigma);
dim   = [T T];
dimp  = dim + 2 * pad;
resadd  = 3;
enlarge = ceil( resadd / 2 );
xvals   = { (1:dimp(1)), (1:dimp(2)) };
% xvals   = { (1:dimp(1))/2, (1:dimp(2))/4 }; % different size of the box in the two directions 
theoryL = LKC_isogauss_theory( FWHM, [ T T ] )

% Mask
mask = true( dim );
mask = logical( pad_vals( mask, pad) );
%%
lat_data = wnfield( mask, nsubj);
LKC_conv_est( lat_data.field, mask, FWHM, resadd)