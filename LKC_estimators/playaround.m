addpath(genpath("/home/drtea/matlabToolboxes/RFTtoolbox/"))

T  = 50;
N  = 30;
rf = noisegen( [T T], N, 5*sqrt(8*log(2)) );
pad = 20;
rf = randn([[50 50]+2*pad, 100]);
mask = ones([T T]);

L = LKCestim_GaussConv( rf, 5*sqrt(8*log(2)), mask, 1, 1, pad );

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

[ Xgrid, Ygrid ] = meshgrid( -siz:dx:siz,  -siz:dx:siz );
x = Xgrid(mask);
y = Ygrid(mask);
DT = delaunayTriangulation( [ x, y ] );
figure(3), clf
triplot(DT,x,y); hold on;
scatter(x,y)

% find boundary
bdy1 = convn(mask, [[0 1];[ -1 0]],'same');
figure(4), clf
imagesc(bdy1);
bdy2 = convn(mask, [[1 1 1];[1 1 1];[1 1 1]],'same')
figure(5), clf
imagesc( bdy2 )

tmp = bdy2 <=3 & bdy2 > 0
tmp 