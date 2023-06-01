%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the LKC_wncfield_theory function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %% D = 1 
%% % Rectangular mask and Kernel
% Create a mask and show it
mask = ones( [ 30 1 ] );
mask = pad_vals( mask, Kernel.truncation );
plot(mask)
title('mask')

%% Stability for increasing resadd
% mask the lattice data
theory_res = 21;
params = ConvFieldParams( repmat(3, [1 1]),...
                          theory_res,...
                          ceil( theory_res / 2 ),...
                          false );

LKC1 = LKC_wncfield_theory( mask, params );

% Show that LKC_conv_est provides similar result (seems to agree)
lat_data = randn( [ size( mask,1 ) 200 ] );
LKChat = LKC_conv_est( lat_data, mask, Kernel, 1, 1 );
LKChat.hatL

% do not mask the lattice data
LKC1 = LKC_wncfield_theory( mask, Kernel, 1, 0 );
LKC3 = LKC_wncfield_theory( mask, Kernel, 3, 0 );
LKC5 = LKC_wncfield_theory( mask, Kernel, 5, 0 );
[ LKC1.L; LKC3.L; LKC5.L ]

% Show that LKC_conv_est provides similar result (seems to agree)
lat_data = randn( [  size( mask,1 ) 200 ] );
LKChat = LKC_conv_est( lat_data, mask, Kernel, 1, 0 );
LKChat.hatL

%% % Complicated Mask and Kernel
% Set Kernel to be isotropic Gaussian with FHWM
Kernel = SepKernel( 1, 20 );

% Create a mask and show it
Sig = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
mask = logical( Sig > 0.02 );
mask = mask(50,:)';
plot( mask ),
title( 'mask' )

%% Stability for increasing resadd
% mask the lattice data
LKC1 = LKC_wncfield_theory( mask, Kernel, 1, 1 );
LKC3 = LKC_wncfield_theory( mask, Kernel, 3, 1 );
LKC5 = LKC_wncfield_theory( mask, Kernel, 5, 1 );
[ LKC1.L; LKC3.L; LKC5.L ]

% Show that LKC_conv_est provides similar result
lat_data = randn( [ size( mask,1 ) 200 ] );
LKChat = LKC_conv_est( lat_data, mask, Kernel, 1, 1 );
LKChat.hatL

% do not mask the lattice data
LKC1 = LKC_wncfield_theory( mask, Kernel, 1, 0 );
LKC3 = LKC_wncfield_theory( mask, Kernel, 3, 0 );
LKC5 = LKC_wncfield_theory( mask, Kernel, 5, 0 );
[ LKC1.L; LKC3.L; LKC5.L ]

% Show that LKC_conv_est provides similar result
lat_data = randn( [ size( mask,1 ) 200 ] );
LKChat = LKC_conv_est( lat_data, mask, Kernel, 1, 0 );
LKChat.hatL

%% %% D = 2 
%% % Rectangular mask and Kernel
% Set Kernel to be isotropic Gaussian with FHWM
Kernel = SepKernel( 2, 3 );

% Create a mask and show it
mask = ones( [ 30 30 ] );
mask = pad_vals( mask, Kernel.truncation );
imagesc(mask)

%% Stability for increasing resadd
% mask the lattice data
LKC1 = LKC_wncfield_theory( mask, Kernel, 1, 1 );
LKC3 = LKC_wncfield_theory( mask, Kernel, 3, 1 );
LKC5 = LKC_wncfield_theory( mask, Kernel, 5, 1 );
[ LKC1.L; LKC3.L; LKC5.L ]

% Show that LKC_conv_est provides similar result (seems to agree)
lat_data = randn( [ size( mask ) 200 ] );
LKChat = LKC_conv_est( lat_data, mask, Kernel, 1, 1 );
LKChat.hatL

% do not mask the lattice data
LKC1 = LKC_wncfield_theory( mask, Kernel, 1, 0 );
LKC3 = LKC_wncfield_theory( mask, Kernel, 3, 0 );
LKC5 = LKC_wncfield_theory( mask, Kernel, 5, 0 );
[ LKC1.L; LKC3.L; LKC5.L ]

% Show that LKC_conv_est provides similar result (seems to agree)
pad_mask = pad_vals( mask, Kernel.truncation );
lat_data = randn( [ size(pad_mask ) 200 ] );
LKChat = LKC_conv_est( lat_data, pad_mask, Kernel, 1, 0 );
LKChat.hatL

%% % Complicated mask and Kernel
% Set Kernel to be isotropic Gaussian with FHWM
Kernel = SepKernel( 2, 20 );

% Create a mask and show it
Sig  = gensig( [1,2], 3, [10,20], [70,100], {[30,20], [50,70]} );
mask = logical( Sig > 0.02 & Sig < 1.1 );
imagesc(mask)

%% Stability for increasing resadd
% mask the lattice data
LKC1 = LKC_wncfield_theory( mask, Kernel, 1, 1 );
LKC3 = LKC_wncfield_theory( mask, Kernel, 3, 1 );
LKC5 = LKC_wncfield_theory( mask, Kernel, 5, 1 );
[ LKC1.L; LKC3.L; LKC5.L ]

% Show that LKC_conv_est provides similar result (seems to agree)
lat_data = randn( [ size( mask ) 200 ] );
LKChat = LKC_conv_est( lat_data, mask, Kernel, 1, 1 );
LKChat.hatL

% do not mask the lattice data
LKC1 = LKC_wncfield_theory( mask, Kernel, 1, 0 );
LKC3 = LKC_wncfield_theory( mask, Kernel, 3, 0 );
LKC5 = LKC_wncfield_theory( mask, Kernel, 5, 0 );
[ LKC1.L; LKC3.L; LKC5.L ]

% Show that LKC_conv_est provides similar result (seems to agree)
lat_data = randn( [ size( mask ) 200 ] );
LKChat = LKC_conv_est( lat_data, mask, Kernel, 1, 0 );
LKChat.hatL

%% %% D = 3
%% % Rectangular mask and Kernel
% Set Kernel to be isotropic Gaussian with FHWM
Kernel = SepKernel( 3, 2 );
FWHM = repmat(3, [1 3]);

% Create a mask
mask = ones( [ 10 10 10 ] );
mask = pad_vals( mask, Kernel.truncation );

%% Stability for increasing resadd
% mask the lattice data
params = ConvFieldParams( FWHM, 1, ceil(1/2), true );
LKC1 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 3, ceil(3/2), true );
LKC3 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 5, ceil(5/2), true );
LKC5 = LKC_wncfield_theory( mask, params );
[ LKC1; LKC3; LKC5 ]

% Show that LKC_conv_est provides similar result (seems to agree)
lat_data = wfield(logical(mask), 50);
params = ConvFieldParams( FWHM, 1, ceil(1/2), true );
LKChat = LKC_latconv_est( lat_data, params );
LKChat

% do not mask the lattice data
params = ConvFieldParams( FWHM, 1, ceil(1/2), false );
LKC1 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 3, ceil(3/2), false );
LKC3 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 5, ceil(5/2), false );
LKC5 = LKC_wncfield_theory( mask, params );
[ LKC1; LKC3; LKC5 ]

% Show that LKC_conv_est provides similar result (seems to agree)
lat_data = wfield(logical(mask), 50);
params = ConvFieldParams( FWHM, 1, ceil(1/2), false );
LKChat = LKC_latconv_est( lat_data, params );
LKChat

%% Mask with holes
x       = linspace(-1,1, 20);
[X,Y,Z] = meshgrid(x,x,x);

mask = exp(-(X.^2 + Y.^2 + Z.^2) / 2 / 5 )>0.911;
%mask(6:8, 10,6:8) = 0;
%mask(12:14,12:14,10) = 0;

% mask the lattice data
params = ConvFieldParams( FWHM, 1, ceil(1/2), true );
LKC1 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 3, ceil(3/2), true );
LKC3 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 5, ceil(5/2), true );
LKC5 = LKC_wncfield_theory( mask, params );
[ LKC1; LKC3; LKC5 ]

% Show that LKC_conv_est provides similar result (seems to agree)
lat_data = wfield(logical(mask), 50);
params = ConvFieldParams( FWHM, 1, ceil(1/2), true );
LKChat = LKC_latconv_est( lat_data, params );
LKChat

% do not mask the lattice data
params = ConvFieldParams( FWHM, 1, ceil(1/2), false );
LKC1 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 3, ceil(3/2), false );
LKC3 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 5, ceil(5/2), false );
LKC5 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 7, ceil(7/2), false );
LKC7 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 9, ceil(9/2), false );
LKC9 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 11, ceil(11/2), false );
LKC11 = LKC_wncfield_theory( mask, params );
[ LKC1; LKC3; LKC5; LKC7; LKC9; LKC11 ]

% Show that LKC_conv_est provides similar result (seems to agree)
lat_data = wfield(logical(mask), 50);
params = ConvFieldParams( FWHM, 1, ceil(1/2), false );
LKChat = LKC_latconv_est( lat_data, params );
LKChat

%% Mask with holes
x = linspace(-1,1, 20);
[X,Y,Z] = meshgrid(x,x,x);

mask = exp(-(X.^2 + Y.^2 + Z.^2) / 2 / 5 )>0.911;
mask(6:8, 10,6:8) = 0;
mask(12:14,12:14,10) = 0;

% mask the lattice data
params = ConvFieldParams( FWHM, 1, ceil(1/2), true );
LKC1 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 3, ceil(3/2), true );
LKC3 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 5, ceil(5/2), true );
LKC5 = LKC_wncfield_theory( mask, params );
[ LKC1; LKC3; LKC5 ]

% Show that LKC_conv_est provides similar result (seems to agree)
lat_data = wfield(logical(mask), 50);
params = ConvFieldParams( FWHM, 1, ceil(1/2), true );
LKChat = LKC_latconv_est( lat_data, params );
LKChat

% do not mask the lattice data
params = ConvFieldParams( FWHM, 1, ceil(1/2), false );
LKC1 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 3, ceil(3/2), false );
LKC3 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 5, ceil(5/2), false );
LKC5 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 7, ceil(7/2), false );
LKC7 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 9, ceil(9/2), false );
LKC9 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 11, ceil(11/2), false );
LKC11 = LKC_wncfield_theory( mask, params );
[ LKC1; LKC3; LKC5; LKC7; LKC9; LKC11 ]

% Show that LKC_conv_est provides similar result (seems to agree)
lat_data = wfield(logical(mask), 50);
params = ConvFieldParams( FWHM, 1, ceil(1/2), false );
LKChat = LKC_latconv_est( lat_data, params );
LKChat

%% % Circular mask and Kernel
% Set Kernel to be isotropic Gaussian with FHWM
Kernel = SepKernel( 3, 2 );
% Create a mask and show it
siz = 3;
dx  = 0.5;
[x,y,z] = meshgrid( -siz:dx:siz, -siz:dx:siz, -siz:dx:siz );
xvals = [x(:), y(:), z(:)]';
h     = reshape( GkerMV( xvals, 5 ), size(x) );
mask  = logical( h > 0.003 );
imagesc( mask(:,:,7) )
clear h

%% Stability for increasing resadd
% mask the lattice data
params = ConvFieldParams( FWHM, 1, ceil(1/2), false );
LKC1 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 3, ceil(3/2), false );
LKC3 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 5, ceil(5/2), false );
LKC5 = LKC_wncfield_theory( mask, params );
params = ConvFieldParams( FWHM, 7, ceil(7/2), false );
LKC7 = LKC_wncfield_theory( mask, params );
[ LKC1; LKC3; LKC5; LKC7 ]

% Show that LKC_conv_est provides similar result (maybe L2 needs a small fix)
lat_data = randn( [ size( mask ) 10 ] );
LKChat = LKC_conv_est( lat_data, mask, Kernel, 1, 1 );
LKChat.hatL

% do not mask the lattice data
LKC1 = LKC_wncfield_theory( mask, Kernel, 1, 0 );
LKC3 = LKC_wncfield_theory( mask, Kernel, 3, 0 );
LKC5 = LKC_wncfield_theory( mask, Kernel, 5, 0 );
[ LKC1.L; LKC3.L; LKC5.L ]

% Show that LKC_conv_est provides similar result (maybe L2 needs a small fix)
lat_data = randn( [ size( mask ) 200 ] );
LKChat = LKC_conv_est( lat_data, mask, Kernel, 1, 0 );
LKChat.hatL