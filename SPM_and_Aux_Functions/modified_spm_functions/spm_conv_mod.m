function [ X, ss ] = spm_conv_mod(X,sx,sy)
% spm_conv is a one or two dimensional convolution of a matrix variable
% in working memory.  It capitalizes on the sparisity structure of the
% problem and the separablity of multidimensional convolution with a Gaussian
% kernel by using one-dimensional convolutions and kernels that are
% restricted to non near-zero values. Need to rewrite your own version of
% this! As this one is not great when the dimensions of the image are low
% in one direction or another.
% 
% If image is mean zero and variance one, scaling the smoothed image by
% 1/sqrt(ss) will return it to unit variance.
%--------------------------------------------------------------------------
% ARGUMENTS
% X     An image which is input as an lx by ly matrix.
% sx 	The FWHM smoothing parameter in the x direction.
% sy  	The FWHM smoothing parameter in the y direction. Note this
%       parameter is optional. If it is omitted then it is deafulted to sx.
%--------------------------------------------------------------------------
% OUTPUT
% X     A matrix representing the image after it has undergone gaussian
%       kernel convolution.
% ss 	The sum of squares of the values that the kernel convolution matrix
%       takes.
%--------------------------------------------------------------------------
% EXAMPLES
% %Smoothing Signal
% Sig = SpheroidSignal([256, 256],20,3,0);
% Img = spm_conv(Sig,10,10);
% surf(Img)
%
% %Smoothing Noise
% Noise = randn([256, 256]); 
% Noise = spm_conv(Noise,5,5);
% surf(Noise)
%--------------------------------------------------------------------------
% SEE ALSO
% spm_smooth
%--------------------------------------------------------------------------
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% Karl Friston
% $Id: spm_conv.m 4195 2011-02-05 18:39:13Z ged $


% assume isometric smoothing
%---------------------------------------------------------------------------
if nargin < 3; sy = sx; end
sx      = abs(sx);
sy      = abs(sy);
[lx, ly] = size(X);

% FWHM -> sigma
%---------------------------------------------------------------------------
sx    = sx/sqrt(8*log(2)) + eps; %Why: so you don't divide by zero.
%Would be better to include an error in this case.
sy    = sy/sqrt(8*log(2)) + eps;

% kernels
%---------------------------------------------------------------------------
Ex    = min([round(6*sx) lx]);  
%It does this as it doesn't bother convolving past a certain point as the
%kernel values will be too low to make a difference.
x     = -Ex:Ex;
kx    = exp(-x.^2/(2*sx^2));
kx    = kx/sum(kx);             %To make each value a proportion of the values around it.
Ey    = min([round(6*sy) ly]);
y     = -Ey:Ey;
ky    = exp(-y.^2/(2*sy^2));
ky    = ky/sum(ky);

if ly > 1 && lx > 1
    [Sx, Sy]  = meshgrid(kx,ky);  %note that meshgrid is the same as ndgrid!
    %Sx: the values the kernel takes in the x direction. Most weight at the
    %given point. Similarly for Sy.
    ss       = sum((Sx(:).*Sy(:)).^2);   %This is the sum of the squares of the values of the kernel here.
elseif ly > 1
    ss = sum(ky.^2);
elseif lx > 1
    ss = sum(kx.^2);
else
    ss = 1;
end

% convolve
%--------------------------------------------------------------------------
%Why does all the flipping go on?
if lx > 1;
    for i = 1:ly
        u      = X(:,i);
        u      = [flipud(u(1:Ex)); u; flipud(u([1:Ex] + lx - Ex))]; 
        %This step takes the first bit and last bit of u and adds them onto
        %the begining and end of u respectively. This is so that the values
        %at the edge are convolved in the right way.
        %Note the original u is a column vector so the new u is too.
        %Why does below use full(u), isn't this unnecessary?
        U      = sparse(conv(full(u),kx)); %Converts to sparse form to save memory space!
        %Check this sparse/full combination, interesting to know when thisI
        %works/doesn't!
        X(:,i) = U([1:lx] + 2*Ex); %Shouldn't this be Ex, not 2*Ex. And similarly for Ey below
        %As you want to grab the middle as this is the relevant bit.
    end
end
if ly > 1; %Ie if there are y dimensions! Ie if X is two dimensional.
    for i = 1:lx
        u      = X(i,:);
        u      = [fliplr(u(1:Ey)) u fliplr(u([1:Ey] + ly - Ey))];
        U      = sparse(conv(full(u),ky));
        X(i,:) = U([1:ly] + 2*Ey);
    end
end
