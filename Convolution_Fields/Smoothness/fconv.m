function [ smoothed_data, ss ] = fconv( data, sep_kern, D, truncation, adjust_kernel )
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
% D = length(Dim); FWHM = 3; D = 3;
% smoothed_spm = zeros(Dim);
% spm_smooth(lat_data, smoothed_spm, FWHM);
% smoothed_fconv = fconv(lat_data, FWHM);
% sigma = FWHM2sigma(FWHM); truncation = ceil(6*sigma);
% smoothed_fconv_spmkern = fconv(lat_data, @(x) spm_smoothkern(FWHM, x), D, truncation );
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
if D > 3
    error('fconv not coded for dimension > 3');
end

if isnumeric(sep_kern)
    FWHM = sep_kern;
    Kernel = @(x) Gker(x, FWHM);
    sigma = FWHM2sigma(FWHM);
    truncation = ceil(4*sigma);
else
    Kernel = sep_kern;
    if nargin < 4
        error('Need to specify truncation')
    end
end

if nargin < 5 || isnan(sum(adjust_kernel(:))) % Allows adjust_kernel to be set as a NaN
    adjust_kernel = zeros(1,D)';
elseif length(adjust_kernel) ~= D
    error('The kernel adjustment must be of the same dimension as the data')
end

% If there are multiple subjects run fconv on each of them
if (D < length(s_data) && D > 1) || (D == 1 && all(s_data > [1,1]))
    smoothed_data = zeros(s_data); nsubj = s_data(end);
    index  = repmat( {':'}, 1, D );
    for J = 1:nsubj
        smoothed_data(index{:}, J) = fconv(squeeze(data(index{:},J)), Kernel, D, truncation, adjust_kernel);
    end
    return
end

% Specify the values at which to evaluate the (1D) kernel at
if length(truncation) > 1 %Allowing the kernel to be evaluated at a user specified set of points in this case
    truncation_vector = truncation;
else % The default: to calculate the kernel at points from -truncation to truncation
    truncation_vector = -truncation:truncation;
end    

% Calculate the kernel at the values in truncation vector (as well as
% including an adjustment if this has been specified)
dthdirectionvector = zeros(D, length(truncation_vector));
kernel_in_direction = zeros(D, length(truncation_vector));
for d = 1:D
    dthdirectionvector(d,:) = fliplr(truncation_vector - adjust_kernel(d));
    kernel_in_direction(d,:) = Kernel(dthdirectionvector(d,:));
end

% Main body, running convolution with a separable kernel in each direction
if D == 1
    smoothed_data = conv(data, kernel_in_direction(1,:),'same');
    ss = sum(dthdirectionvector(1,:).^2); % Calculates the sum of the squares of the kernel
elseif D == 2
    xside = kernel_in_direction(1,:)'; %The kernel is the x direction (x corresponds to the rows of the matrix)
    yside = kernel_in_direction(2,:); % The kernel in the y direction (y corresponds to the rows of the matrix)
    smoothed_data = convn(data, xside, 'same'); %Smooth in the x direction
    smoothed_data = convn(smoothed_data, yside, 'same'); %Smooth in the y direction
%     smoothed_data = convn(udside, smoothed_data, 'same');  %Smooth in the up down direction
    
    [sx,sy] = meshgrid(dthdirectionvector(1,:),dthdirectionvector(2,:)); % Calulate the kernel values everywhere
    ss = sum((sx(:).*sy(:)).^2); % This is the sum of the squares of the kernel in 
                            % the truncation box, since the kernel is by assumption separable.
elseif D == 3
    xside = kernel_in_direction(1,:)'; %The kernel is the x direction
    yside = kernel_in_direction(2,:); %The kernel is the y direction
    zside = zeros(1,1,length(kernel_in_direction(3,:)));  
    zside(1,1,:) = kernel_in_direction(3,:); %The kernel is the z direction
    smoothed_data = convn(data, xside, 'same'); %Smooth in the x direction
    smoothed_data = convn(smoothed_data, yside, 'same'); %Smooth in the y direction
    smoothed_data = convn(smoothed_data, zside, 'same'); %Smooth in the z direction
     
    [sx,sy,sz] = meshgrid(dthdirectionvector(1,:),dthdirectionvector(2,:),dthdirectionvector(3,:));
    ss = sum((sx(:).*sy(:).*sz(:)).^2); % This is the sum of the squares of the kernel in 
                            % the truncation box, since the kernel is by assumption separable.
else
    error('fconv not coded for dimension > 3')
end

end

% Deprecated
% if D == 1
%     smoothed_data = conv(smoothed_data, kernel_in_direction(d,:),'same'); 
%     % Could use convn here, but this way gives the same result no matter 
%     % whether data is a row or a column vector
% else
%     for d = 1:D
%         smoothed_data = convn(smoothed_data, kernel_in_direction(d,:),'same');
%     end
% end
