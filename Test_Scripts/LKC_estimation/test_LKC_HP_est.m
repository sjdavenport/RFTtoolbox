%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%    This script tests the LKC_HP_est.m function
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
dim   = 100;
mask  = ones( [dim 1] );
FWHM  = 12;

% generate random noise
Y = noisegen( dim, 10, FWHM, 0 );

% compute HPE for one subject. Note that normalize needs to be 0! It is
% automatically correctly chosen if not specified. All values give the same
% result!
HPE_C = LKC_HP_est( Y(:,1), 1, 1, "C" );
HPE_m = LKC_HP_est( Y(:,1), 1, 1, "matlab" );
HPE_a = LKC_HP_est( Y(:,1), mask );

if ~( HPE_C.hatL == HPE_m.hatL && HPE_C.hatL == HPE_a.hatL...
        && HPE_a.hatL == HPE_m.hatL )
    error( "The values of the 3 estimates are unequal!" )
end
 
% compute HPE for several subject. Here we can choose to normalize, which
% is recommended.
for normalize = [0 1]
    HPE_C = LKC_HP_est( Y, mask, 1, 1, "C" );
    HPE_m = LKC_HP_est( Y, mask, 1, 1, "matlab" );
    HPE_a = LKC_HP_est( Y, mask ); % same as HPE_m using "C"

    if ~( HPE_C.hatL == HPE_m.hatL && HPE_C.hatL == HPE_a.hatL...
            && HPE_a.hatL == HPE_m.hatL )
        error( "The values of the 3 estimates are unequal!" )
    end
end

% compute HPE for several subject. Here normalize does not have any
% since it is required to normalize and done automatically. Note that bHPE
% requires several subjects
rng(1)
bHPE_C = LKC_HP_est( Y, mask, 1e3, 0, "C" );
rng(1)
bHPE_m = LKC_HP_est( Y, mask, 1e3, 0, "matlab" );

if bHPE_C.hatL ~= bHPE_m.hatL
    error( "The values of the estimates are unequal!" )
end

rng(1)
bHPE_C = LKC_HP_est( Y, mask, 1e3, 1, "C" );
rng(1)
bHPE_m = LKC_HP_est( Y, mask, 1e3, 1, "matlab" );

if bHPE_C.hatL ~= bHPE_m.hatL
    error( "The values of the estimates are unequal!" )
end


%% %-----------------------------------------------------------------------
%    Test section D = 2
%--------------------------------------------------------------------------
dim   = [ 50 50 ];
mask  = ones( dim );
FWHM  = sigma2FWHM( 5 );
mask( 25, 25 ) = -Inf;

% generate random noise
Y = noisegen( dim, 10, FWHM, 0 );

% compute HPE for one subject. Note that normalize needs to be 0! It is
% automatically correctly chosen if not specified. All values give the same
% result!
HPE_C = LKC_HP_est( Y(:,:,1), mask, 1, 0, "C" );
HPE_m = LKC_HP_est( Y(:,:,1), mask, 1, 0, "matlab" );
HPE_a = LKC_HP_est( Y(:,:,1), mask );

if ~( all( HPE_C.hatL == HPE_m.hatL ) && all( HPE_C.hatL == HPE_a.hatL )...
        && all( HPE_a.hatL == HPE_m.hatL ) )
    error( "The values of the 3 estimates are unequal!" )
end
 
% compute HPE for several subject. Here we can choose to normalize, which
% is recommended.
for normalize = [0 1]
    HPE_C = LKC_HP_est( Y, mask, 1, 1, "C" );
    HPE_m = LKC_HP_est( Y, mask, 1, 1, "matlab" );
    HPE_a = LKC_HP_est( Y, mask ); % same as HPE_m using "C"

    if ~( all( HPE_C.hatL == HPE_m.hatL ) && all( HPE_C.hatL == HPE_a.hatL )...
            && all( HPE_a.hatL == HPE_m.hatL ) )
        error( "The values of the 3 estimates are unequal!" )
    end
end

% compute HPE for several subject. Here normalize does not have any
% since it is required to normalize and done automatically. Note that bHPE
% requires several subjects
for normalize = [0 1]
    rng(1)
    bHPE_C = LKC_HP_est( Y, mask, 1e3, 0, "C" );
    rng(1)
    bHPE_m = LKC_HP_est( Y, mask, 1e3, 0, "matlab" );

    if any( bHPE_C.hatL ~= bHPE_m.hatL )
        error( "The values of the estimates are unequal!" )
    end
end


%% %-----------------------------------------------------------------------
%    Test section D = 3
%--------------------------------------------------------------------------
dim   = [ 35 35 35 ];
mask  = ones( dim );
FWHM  = 5;
mask( 25, 25 ) = 1;

% generate random noise
Y = noisegen( dim, 10, FWHM );

% compute HPE for one subject. Note that normalize needs to be 0! It is
% automatically correctly chosen if not specified. All values give the same
% result!
HPE_C = LKC_HP_est( Y(:,:,:,1), mask, 1, 0, "C" );
HPE_m = LKC_HP_est( Y(:,:,:,1), mask, 1, 0, "matlab" );
HPE_a = LKC_HP_est( Y(:,:,:,1), mask );

% compute HPE for several subject. Here we can choose to normalize, which
% is recommended.
HPE_C = LKC_HP_est( Y, mask, 1, 0, "C" );
HPE_m = LKC_HP_est( Y, mask, 1, 1, "matlab" );
HPE_a = LKC_HP_est( Y, mask ); % same as HPE_m using "C"

% compute HPE for several subject. Here normalize does not have any
% since it is required to normalize and done automatically. Note that bHPE
% requires several subjects. Note value is unexpectedly low. Need to figure
% out why. Never use "matlab" for bHPE since it is extremely slow.
bHPE_C = LKC_HP_est( Y, mask, 3e3, 0, "C" );
