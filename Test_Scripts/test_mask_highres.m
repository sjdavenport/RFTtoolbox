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

mask_hr = mask_highres( mask, resAdd, enlarge );

%% 2D
% resolution added
resAdd  = 3;
% create a mask and show it
Sig  = gensig([1,2], 3, [10,20], [100,150], {[40,30], [70,120]});
mask = logical( Sig > 0.02 & Sig < 1.1 );
imagesc( mask )
clear Sig

% enlarged domain
mask_hr = mask_highres( mask, resAdd, 1, show_plot );
% non enlarged domain
mask_hr = mask_highres( mask, resAdd, 0, show_plot );

% resolution added
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
