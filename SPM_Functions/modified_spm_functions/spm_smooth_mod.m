function [ss] = spm_smooth(P,Q,s,dtype)
% 3 dimensional convolution of an image
% FORMAT spm_smooth(P,Q,s,dtype)
% P     - image(s) to be smoothed (or 3D array)
% Q     - filename for smoothed image (or 3D array)
% s     - [sx sy sz] Gaussian filter width {FWHM} in mm (or edges)
% dtype - datatype [Default: 0 == same datatype as P]
% 
% ss     - Sums of squares of kernal
%__________________________________________________________________________
%
% spm_smooth is used to smooth or convolve images in a file (maybe).
%
% The sum of kernel coeficients are set to unity.  Boundary
% conditions assume data does not exist outside the image in z (i.e.
% the kernel is truncated in z at the boundaries of the image space). s
% can be a vector of 3 FWHM values that specifiy an anisotropic
% smoothing.  If s is a scalar isotropic smoothing is implemented.
%
% If Q is not a string, it is used as the destination of the smoothed
% image.  It must already be defined with the same number of elements
% as the image.
%
% If image is mean zero and variance one, scaling smoothed image by
% 1/sqrt(ss) will return it to unit variance.
%
%__________________________________________________________________________
% Copyright (C) 1994-2011 Wellcome Trust Centre for Neuroimaging
% John Ashburner & Tom Nichols
% $Id: spm_smooth.m 4419 2011-08-03 18:42:35Z guillaume $
% T. Nichols mode: Returns the sum of squares of the kernel
%
% SEE ALSO
% spm_vol, smooth1, spm_conv_vol, spm_create_vol, spm_smoothkern

%-Parameters & Arguments
%--------------------------------------------------------------------------
if numel(s) == 1, s = [s s s]; end
if nargin < 4, dtype = 0; end

if ischar(P), P = spm_vol(P); end

if isstruct(P) %This allows for multiple Ps to be read in!
    for i= 1:numel(P)
        ss1=smooth1(P(i),Q,s,dtype);
    end
else
    %ss1=smooth1(P,Q,s,dtype);
end
% if nargout>0
%   ss = ss1;
% end

% if isstruct(P)
%     VOX = sqrt(sum(P.mat(1:3,1:3).^2));
% else
%     VOX = [1 1 1];
% end
% 
% %-Compute parameters for spm_conv_vol
% %--------------------------------------------------------------------------
% s  = s./VOX;                        % voxel anisotropy
 s1 = s/sqrt(8*log(2));              % FWHM -> Gaussian parameter

x  = round(6*s1(1)); x = -x:x; x = spm_smoothkern(s(1),x,1); x  = x/sum(x);
y  = round(6*s1(2)); y = -y:y; y = spm_smoothkern(s(2),y,1); y  = y/sum(y);
z  = round(6*s1(3)); z = -z:z; z = spm_smoothkern(s(3),z,1); z  = z/sum(z);

i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;

if isstruct(Q), Q = spm_create_vol(Q); end

spm_conv_vol(P,Q,x,y,z,-[i,j,k]);

[sx,sy,sz] = meshgrid(x,y,z);
ss         = sum((sx(:).*sy(:).*sz(:)).^2);
