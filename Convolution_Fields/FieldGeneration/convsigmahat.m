function [ sigma ] = convsigmahat( tval, lat_data, Kernel, mask, truncation, xvals_vecs, meanfn )
% Returns just the sigma value from applying applyconvfield_t, see applyconvfield_t
% for information on the rest of the functions.
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

if ~exist('truncation', 'var')
    truncation = -1;
end
if isnumeric(Kernel)
    if truncation == -1
        sigma = FWHM2sigma(Kernel);
        truncation = round(4*sigma); % Default truncation
    end
    Kernel = @(x) GkerMV(x,Kernel); % The multivariate Gaussian kernel
else
    %     Need to work out what to do here for general kernels
    if truncation == -1
        % The -1 selection is only set up for isotropic Gaussian kernels atm
        truncation = 0; 
    end
end

% Default is to take the mean to be 0
if ~exist('meanfn', 'var')
    meanfn = @(x) 0;
end

Ldim = size(lat_data);
D = length(Ldim) - 1;

if ~exist('xvals_vecs', 'var')
    xvals_vecs = {1:Ldim(1)}; %The other dimensions are taken case of below.
end

if ~iscell(xvals_vecs)
    xvals_vecs = {xvals_vecs};
end

if length(xvals_vecs) < D
    increm = xvals_vecs{1}(2) - xvals_vecs{1}(1);
    for d = (length(xvals_vecs)+1):D
        xvals_vecs{d} = xvals_vecs{1}(1):increm:(xvals_vecs{1}(1)+increm*(Ldim(d)-1));
    end
end

if ~exist('mask', 'var') %Default mask
    if D == 1
        mask = ones([1,Ldim(1)]);
%         mask = ones([1,Ldim(1:D)]);
    else
        mask = ones(Ldim(1:D));
    end
else 
    if D == 1
       % Need to ensure that mask is a row vector for use in applyconvfield
       if size(mask, 2) == 1
           mask = mask'; 
       end
    end
end

nsubj = size(lat_data,D+1);
nevals = size(tval,2);
Xcfields_at_tval = zeros(nevals, nsubj);

% if size(data,1) ~= length(xvalues_at_voxels)
%     error('The dimensions of data and xvalues_at_voxels do not match up')
% end

if D == 1
    for I = 1:nsubj
        Xcfields_at_tval(:, I) = applyconvfield(tval, lat_data(:,I)', Kernel, mask, truncation, xvals_vecs ) + meanfn(tval);
    end
elseif D == 2
    for I = 1:nsubj
        Xcfields_at_tval(:, I) = applyconvfield(tval, squeeze(lat_data(:,:,I)), Kernel, mask, truncation, xvals_vecs )  + meanfn(tval);
    end
elseif D == 3
    for I = 1:nsubj
        Xcfields_at_tval(:, I) = applyconvfield(tval, squeeze(lat_data(:,:,:,I)), Kernel, mask, truncation, xvals_vecs )  + meanfn(tval);
    end
else
    error('Not coded for D > 3')
end

[~,~,sigma] = mvtstat(Xcfields_at_tval, 0);

sigma = sigma';

end

