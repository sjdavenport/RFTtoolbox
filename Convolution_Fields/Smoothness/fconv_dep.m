function [ smoothed_data, ss ] = fconv_dep( data, sep_kern, D, adjust_kernel, truncation )
% FCONV( data, sep_kern, truncation ) provides a faster implementation for
% smoothing data using a separable kernel (e.g. an isotropic Gaussian
% kernel).
%--------------------------------------------------------------------------
% ARGUMENTS
% data          a Dim by nsubj array of data
% sep_kern      a function handle giving a separable kernel. If this is 
%               instead numeric fconv smoothes with an isotropic Gaussian 
%               kernel with sep_kern as the FWHM (see EXAMPLES section)
% D             the dimension
% truncation    the truncation of the Kernel to use (if using a the
%--------------------------------------------------------------------------
% OUTPUT
% smoothed_data     the smoothed data
% ss                the sum of squares of the kernel (useful for ensuring
%                   variance 1 isotropic fields
%--------------------------------------------------------------------------
% EXAMPLES
% %1D
% lat_data = normrnd(0,1,1,100); FWHM = 3;
% smoothed_fconv = fconv(lat_data, FWHM);
% smoothed_spm = spm_conv(lat_data,FWHM);
% plot(smoothed_spm); hold on; plot(smoothed_fconv)
% legend('spm\_conv', 'fconv') 
%
% %1D multiple subjects
% nvox = 100; nsubj = 2; lat_data = normrnd(0,1,nvox,nsubj); FWHM = 3; D = 1;
% smoothed_fconv = fconv(lat_data, FWHM, D)
% smoothed_spm = zeros(nvox, nsubj);
% for n = 1:nsubj
%     smoothed_spm(:,n) = spm_conv(lat_data(:,n),FWHM);
% end
% plot(smoothed_spm, 'color',[0.85 0.325 0.0980]); hold on;
% plot(smoothed_fconv, 'color',[0 0.447 0.7410]);
%
% %%2D temp
% lat_data = normrnd(0,1,25,25); FWHM = 3;
% smoothed_fconv = fconv(lat_data, FWHM); 
% smoothed_fconvdep = fconv_dep(lat_data,FWHM)
% subplot(2,1,1)
% surf(smoothed_fconv)
% title('fconv')
% subplot(2,1,2)
% surf(smoothed_fconvdep)
% title('fconv dep')
%
% % 2D 
% lat_data = normrnd(0,1,25,25); FWHM = 3;
% smoothed_fconv = fconv(lat_data, FWHM); 
% smoothed_spm = spm_conv(lat_data,FWHM)
% subplot(2,1,1)
% surf(smoothed_fconv)
% title('fconv')
% subplot(2,1,2)
% surf(smoothed_fconv)
% title('SPM\_conv')
%
% % 3D
% Dim = [50,50,50]; lat_data = normrnd(0,1,Dim); halfDim = Dim(1)/2;
% D = length(Dim);
% smoothed_spm = zeros(Dim);
% spm_smooth(lat_data, smoothed_spm, FWHM);
% smoothed_fconv = fconv(lat_data, FWHM);
% sigma = FWHM2sigma(FWHM); truncation = ceil(6*sigma);
% smoothed_fconv_spmkern = fconv(lat_data, @(x) spm_smoothkern(FWHM, x), 3, truncation );
% smoothed_cfield = convfield( lat_data, FWHM, 1, D, 0, 0);
% plot(1:Dim(1),smoothed_fconv(:,halfDim,halfDim))
% hold on 
% plot(1:Dim(1),smoothed_spm(:,halfDim,halfDim))
% plot(1:Dim(1),smoothed_cfield(:,halfDim,halfDim), '--')
% plot(1:Dim(1),smoothed_fconv_spmkern(:,halfDim,halfDim), '--')
% legend('fconv', 'SPM', 'convfield', 'fconv_smoothkern')
%
% plot(-truncation:truncation, spm_smoothkern(FWHM, -truncation:truncation))
% hold on
% plot(-truncation:truncation, GkerMV(-truncation:truncation, FWHM))
%
% % Compare speed to spm_smooth (much faster)
% Dim = [50,50,50]; lat_data = normrnd(0,1,Dim);
% tic; fconv(lat_data, FWHM); toc
% tic; smoothed_spm = zeros(Dim);
% tt = spm_smooth_mod(lat_data, smoothed_spm, FWHM); toc
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
s_data = size(data);
if nargin < 3
    s_data = size(data);
    D = length(s_data);
    if D == 2
        if s_data(1) == 1 
            D = 1;
        elseif s_data(2) == 1
            D = 1;
            data = data';
        end
    end
end
if isnumeric(sep_kern)
    FWHM = sep_kern;
    Kernel = @(x) GkerMV(x, FWHM);
    sigma = FWHM2sigma(FWHM);
    truncation = ceil(4*sigma);
else
    Kernel = sep_kern;
    if nargin < 5
        error('Need to specify truncation')
    end
end

if nargin < 4 || isnan(adjust_kernel)
    adjust_kernel = zeros(1,D)';
elseif length(adjust_kernel) ~= D
    error('The kernel adjustment must be of the same dimension as the data')
end

if nargin > 3
    if length(truncation) > 1 %Allowing the kernel to be evaluated at a user specified set of points in this case
        gridside = Kernel(truncation);
    else % The default: to calculate the kernel at points from -truncation to truncation
        gridside = Kernel(-truncation:truncation);
    end
else
    gridside = Kernel(-truncation:truncation);
end

% If there are multiple subjects run fconv on each of them
if (D < length(s_data) && D > 1) || (D == 1 && all(s_data > [1,1]))
    smoothed_data = zeros(s_data); nsubj = s_data(end);
    index  = repmat( {':'}, 1, D );
    for J = 1:nsubj
        smoothed_data(index{:}, J) = fconv_dep(squeeze(data(index{:},J)), Kernel, D, adjust_kernel, truncation);
    end
    return
end

% Main body, running convolution with a separable kernel in D = 1,2,3
if D == 1
    smoothed_data = conv(data, gridside,'same');
    ss = sum(gridside.^2);
elseif D == 2
    udside = gridside';
    smoothed_data = convn(data, gridside, 'same');
    smoothed_data = convn(smoothed_data, udside, 'same');
    
    [sx,sy] = meshgrid(gridside,gridside);
    ss = sum((sx(:).*sy(:)).^2);
elseif D == 3
    udside = gridside';
    inside = zeros(1,1,length(gridside));
    inside(1,1,:) = gridside;
    smoothed_data = convn(data, gridside, 'same');
    smoothed_data = convn(smoothed_data, udside, 'same');
    smoothed_data = convn(smoothed_data, inside, 'same');
    
    [sx,sy,sz] = meshgrid(gridside,gridside,gridside);
    ss = sum((sx(:).*sy(:).*sz(:)).^2);
else
    error('fconv not coded for dimension > 3')
end

end

