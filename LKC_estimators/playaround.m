clear all
close all

addpath(genpath("/home/drtea/matlabToolboxes/RFTtoolbox/"))

T  = 50;
N  = 30;
D  = 2;

dim  = repmat(T, [1 D]);
pad = 20;
rf = randn([dim+2*pad, 100]);
mask = ones(dim);

L1 = LKCestim_GaussConv( rf, 5*sqrt(8*log(2)), mask, 1, pad );
L1

mask2 = zeros(dim+2*pad);
mask2( (pad+1):(end-pad), (pad+1):(end-pad) ) = mask;
rf = noisegen( 2*dim+2*pad, N, 5*sqrt(8*log(2)) );

L2 = LKCestim_HPE( rf, D, mask2, 1, "C" );
L2.hatn

%



%%
siz = 5;
FWHM = 3;
dx = 1;

siz = ceil( 4*FWHM2sigma(FWHM) );
[x,y] = meshgrid( -siz:dx:siz, -siz:dx:siz );
xvals = [x(:), y(:)]';

% convolution kernels to be used ith convn
mask   = reshape( GkerMV( xvals, FWHM ), size(x) );
figure(1)
imagesc(mask), colorbar;
mask = boolean(mask>2e-3);
figure(2)
imagesc(mask)
figure(3)
imagesc(dilateSet( ~mask ))

bdry_mask = dilateSet( mask ) | dilateSet(~mask)
%%

mask(:,1:end-1)

[ Xgrid, Ygrid ] = meshgrid( -siz:dx:siz,  -siz:dx:siz );
x = Xgrid(mask);
y = Ygrid(mask);
DT = delaunayTriangulation( [ x, y ] );
figure(3), clf
triplot(DT,x,y); hold on;
scatter(x,y)

% find boundary
bdy1 = convn(mask, [[1 1 1];[ 1 1 1];[ 1 1 1]],'same');
figure(4), clf
imagesc(bdy1);
bdy2 = convn(mask, [[1 1 1];[1 1 1];[1 1 1]],'same')
figure(5), clf
imagesc( bdy2 )

tmp = bdy2 <=3 & bdy2 > 0
tmp 



%% %%%%%%%%%%%% memory usage
clear all
close all
domDim = [50 50];
D      = length(domDim);
N      = 5000;

FWHM   = 5*2*sqrt(2*log(2));
pad    = round(FWHM/(2*sqrt(2*log(2)))*5);

Y      = randn([domDim+2*pad, N]);

resAdd   = 2;
remove2 = pad * ( 1 + resAdd);

dx       = 1/(resAdd+1);
domDimhr = ( domDim+2*pad - 1 ) * resAdd + domDim+2*pad;
switch D
    case 3
        % increase the resolution of the raw data by introducing zeros
        Y2 = zeros( [ domDimhr, N ] );
        Y2( 1:( resAdd + 1 ):end, 1:( resAdd + 1 ):end, 1:( resAdd + 1 ):end, : ) = Y;
        % grid for convolution kernel
        siz = ceil( 4*FWHM2sigma(FWHM) );
        [x,y,z] = meshgrid( -siz:dx:siz, -siz:dx:siz, -siz:dx:siz );
        xvals = [x(:), y(:), z(:)]';

        % convolution kernels to be used ith convn
        h   = reshape( GkerMV( xvals, FWHM ), size(x) );
    case 2
        % increase the resolution of the raw data by introducing zeros
        Y2 = zeros( [ domDimhr, N ] );
        Y2( 1:( resAdd + 1 ):end, 1:( resAdd + 1 ):end, : ) = Y;
        % grid for convolution kernel
        siz = ceil( 4*FWHM2sigma(FWHM) );
        [x,y] = meshgrid( -siz:dx:siz, -siz:dx:siz );
        xvals = [x(:), y(:)]';
        
        % convolution kernels to be used ith convn
        h   = reshape( GkerMV( xvals, FWHM ), size(x) );
end

h = h / sqrt(sum(h(:).^2));

% get the convolutional field
smY  = convn( Y2, h, 'same' );

smY = smY( (remove2 + 1):(end-remove2), (remove2 + 1):(end-remove2), : );

figure(1)
imagesc( sqrt(var(smY,0,3))), colorbar
figure(2)
imagesc(smY(1:(resAdd+1):end,1:(resAdd+1):end,4)), colorbar

% D=3: res = 1, 0.2sec/field, res = 3, 1.6sec/field, dim=[50 50 50]
sm = size(smY)
tic
L2 = LKCestim_HPE( smY, D, ones(sm(1:D)), 5e2, "C" );
toc
L2.hatn

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% play around with 3D LKCconv estimator

clear all
close all

addpath(genpath("/home/drtea/matlabToolboxes/RFTtoolbox/"))

T  = 30;
N  = 30;
D  = 3;

FWHM = 3

dim  = repmat(T, [1 D]);
pad = ceil(1.7*FWHM);
rf = randn([dim+2*pad, N]);
mask = ones(dim);

tic
L1 = LKCestim_GaussConv( rf, FWHM, mask, 0, pad );
toc
L1

% Theory
alpha = 1/(4*(FWHM/sqrt(8*log(2)))^2);
LKC = [3*(T-1); 3*(T-1)^2; (T-1)^3] .* [sqrt(2*alpha); 2*alpha; (2*alpha)^(3/2)] ;

LKC(3)

%% understanding Sam's covering function
clear all
close all

addpath( genpath( "/home/drtea/matlabToolboxes/RFTtoolbox/" ) )
dim = 30*ones(1, 2);
data = randn( [ dim 5000 ] );
mask = ones(dim);

mean_record_coverage( data, 50, 3, mask )
