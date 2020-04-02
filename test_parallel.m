clear all
close all

addpath(genpath("/home/drtea/matlabToolboxes/RFTtoolbox"))
addpath(genpath("/home/drtea/matlabToolboxes/spm12"))

% generate 3D field 
FWHM = 5*2*sqrt(2*log(2));
nsubj = 75;
Dim = [50, 50];
D = length(Dim);
f = noisegen( Dim, nsubj, FWHM, 0);
mask = ones(Dim);

% number of bootstraps used
Mboot   = 1e3;
locWork = 2; 

delete(gcp) 
parpool( locWork );
% slow parallel version
rng(1)
tic;
 LKC1 = LKCestim_HPE( f, D, mask, Mboot, 0, num2str(locWork) );
toc;
% quick parallel version
rng(1)
tic;
 LKC2 = LKCestim_HPE2( f, D, mask, Mboot, 0, num2str(locWork) );
toc;
delete(gcp)
% check values
mean( (LKC1.hat1(:)-LKC2.hat1(:)).^2 )
 
% normal version 1
rng(1)
tic;
 LKC3 = LKCestim_HPE( f, D, mask, Mboot );
toc;
% check values
mean( (LKC2.hat1(:)-LKC3.hat1(:)).^2 )

% normal version 2
rng(1)
tic;
LKC4 = LKCestim_HPE( f, D, mask, Mboot, 0, "C" );
toc;
% check values
mean( (LKC2.hat1(:)-LKC4.hat1(:)).^2 )
