%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%    This script tests the mask_highres function
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% prepare workspace
clear all
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1D
% resolution added
resAdd  = 1;
enlarge = 0;
% create a mask and show it
Sig  = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
mask = logical( Sig > 0.02 );
mask = mask(50,:)';
plot( mask )
clear Sig

mask_hr = mask_highres( mask, resAdd, enlarge, 1 ); % What is shown exactly in these plots??

%% 2D issue (there is a bug here!)
Dim = [3,3]; mask = true(Dim); mask(2,1) = 0; resadd = 2;
mask_hr = mask_highres(mask, resadd, 0);
imagesc(mask_hr)

mask_hr = mask_highres(mask, resadd, ceil(resadd/2),1);
imagesc(mask_hr)

% This isn't right, in the enlarged version the voxels should be the same
% size as in the original version. They shouldn't change shape or size?
% Voxels should still be square but clearly not here!
% The "enlarged" version should match the original mask not reshape it!

%% Dividing things into 8
mask = true(2);
mask_highres(mask, 2, 1)

mask = true(3); mask(2,2) = 0;
mask_highres(mask, 3, 0)


%% 2D
% resolution added
resAdd  = 3;
% create a mask and show it
Sig  = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
mask = logical( Sig > 0.02 & Sig < 1.1 );
imagesc( mask )
clear Sig

% enlarged domain
mask_hr = mask_highres( mask, resAdd, ceil(resadd/2),1); %This is definitely not what we want for the enlarged mask?
imagesc(mask_hr)
%% non enlarged domain
mask_hr = mask_highres( mask, resAdd, 0 );
imagesc(mask_hr)

%% resolution added
resAdd  = 3;
enlarge = 1;
% create a mask and show it
mask = logical( ones(50,50) );
mask(:,50) = 0;
imagesc( mask )
clear Sig

% enlarged domain
mask_hr = mask_highres( mask, resAdd, 1, show_plot );
% non enlarged domain
mask_hr = mask_highres( mask, resAdd, 0, show_plot );


%% 3D
% resolution added
resAdd  = 3;
% create a mask and show it
siz = 13;
dx  = 0.25;
[x,y,z] = meshgrid( -siz:dx:siz, -siz:dx:siz, -siz:dx:siz );
xvals = [x(:), y(:), z(:)]';
h     = reshape( GkerMV( xvals, 5 ), size(x) );
mask  = logical( h > 0.002 );
imagesc( mask(:,:,53) )
clear h

% enlarged domain
mask_hr = mask_highres( mask, resAdd, 1, 0 );
% non enlarged domain
mask_hr = mask_highres( mask, resAdd, 0, 0 );
