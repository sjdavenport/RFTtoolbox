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
Mboot = 5e3;


delete(gcp) 
parpool( 2 );
tic;
 LKC = LKCestim_HPE( f, D, mask, Mboot, 0, "2" );
toc;
L1 = LKC.hatn;
 
delete(gcp('nocreate')) 
parpool( 2 );
tic;
 LKC = LKCestim_HPE2( f, D, mask, Mboot, 0, "2" );
toc;
L3 = LKC.hatn;
 
 delete(gcp)
tic;
 LKC = LKCestim_HPE( f, D, mask, Mboot, 0, "C" );
toc;
 L2 = LKC.hatn;

 % difference between estimates obtained by parallelization or not
dL = L1-L2

parTimeE - parTimeB
norTimeE - norTimeB