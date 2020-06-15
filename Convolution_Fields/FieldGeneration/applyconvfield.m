function [field_vals, ss] = ...
       applyconvfield(tval, lat_data, Kernel, mask, truncation, xvals_vecs)
% APPLYCONVFIELD(tval, lat_data, Kernel, mask, truncation, xvals_vecs)
% calculates the convolution field at points specified by tval
%--------------------------------------------------------------------------
% ARGUMENTS
% tval      the t values (an ndim = D by nvalues matrix) at which to evaluate 
%           the field.
% Y         the lattice field an array of size Dim.
% Kernel    the smoothing Kernel (as a function). If Kernel is a postive
%           number then an isotropic Gaussian kernel with dimension D and
%           FWHM = Kernel is used
% truncation    a window around the points at which to evaluate the kernel
%               setting this parameter allows for quicker computation of
%           the convolution field and has very litle effect on the values
%           of the field for kernels that have light tails such as the
%           Gaussian kernel. Default (which is recorded by setting
%           truncation = -1) results in a truncation of 4*FWHM2sigma(Kernel);
%           Setting truncation = 0 results in no truncation at all.
% xvals_vecs    a D-dimensional cell array whose entries are vectors giving the
%               xvalues at each each dimension along the lattice. It assumes
%               a regular, rectangular lattice (though within a given
%               dimension the voxels can be spaced irregularly).
%               I.e suppose that your initial lattice grid is a
%               4by5 2D grid with 4 voxels in the x direction and 5 in
%               the y direction. And that the x-values take the values:
%               [1,2,3,4] and the y-values take the values: [0,2,4,6,8].
%               Then you would take xvals_vecs = {[1,2,3,4], [0,2,4,6,8]}.
%               The default is to assume that the spacing between the
%               voxels is 1. If only one xval_vec direction is set the
%               others are taken to range up from 1 with increment given by
%               the set direction.
% mask      a 0/1 array with the same number of dimensions as Y that serves
%           as a mask for the data. The default (if no mask is included) is
%           not to mask the data, i.e. setting mask = ones(Dim);
%--------------------------------------------------------------------------
% OUTPUT
% field_vals    an nvalues length vector whose entries are the field 
%               evaluated at the locations specified by xvals
%--------------------------------------------------------------------------
% EXAMPLES  
% %% %% 1D examples
% %% Simple 1D Example
% Y = [1,2,3,4];
% tval = 2; FWHM = 3;
% applyconvfield(tval, Y, FWHM)
% 
% %% 1D Comparison with convfield
% nvox = 50; lat_data = randn(1, nvox); FWHM = 3; resadd = 2;
% acf = @(tval) applyconvfield(tval, lat_data, FWHM);
% [cfield, xvals_vecs] = convfield( lat_data, FWHM, resadd, 1);
% 
% plot(xvals_vecs{1}, acf(xvals_vecs{1}))
% hold on
% plot(xvals_vecs{1}, cfield, '--');
% 
% %% 1D Comparison with spm\_conv (if installed)
% nvox = 100; lat_field = normrnd(0,1,nvox,nsubj);
% field_at_voxels = applyconvfield(1:nvox, mean(lat_field,2)', FWHM);
% [~, ss] = applyconvfield(nvox/2, lat_field(:,1)', FWHM);
% plot(field_at_voxels/sqrt(ss))
% hold on
% [smoothfield, ss_spm] = spm_conv_mod(mean(lat_field,2), FWHM);
% plot(smoothfield'/sqrt(ss_spm));
% 
% %% 1D masking
% FWHM = 4; nvox = 50; D = 1;
% mask = true([1,nvox]); mask(20:40) = 0;
% lat_data = randn(1,nvox); resadd = 2;
% field = @(tval) applyconvfield( tval, lat_data, FWHM, mask );
% mfield = @(x) mask_field(x, mask);
% masked_field = @(x) mfield(x).*field(x);
% 
% [cfield, xvals_vecs] = convfield( lat_data, FWHM, resadd, D );
% 
% plot(xvals_vecs{1}, masked_field(xvals_vecs{1}))
% hold on
% plot(xvals_vecs{1}, cfield, '--');
% legend('masked field', 'unmasked field')
% 
% %% %% 2D examples
% %% Simple 2D example
% Y = [1,2;3,4];
% tval = [1.5,2,4; 3.4,2, 5]; FWHM = 3;
% applyconvfield(tval, Y, FWHM)
% 
% %% Demonstrating the need for truncation (and for using applyconvfield)
% % Need to truncate for speed, else  things are really slow!
% FWHM = 3; lat_data = normrnd(0,1,1000,1000);
% tic; applyconvfield([500,500]', lat_data, FWHM); toc
% tic; applyconvfield([500,500]', lat_data, FWHM, 10); toc
% tic; convfield(lat_data, FWHM, 0, 2); toc
% 
% %% %% 3D examples
% %% Truncation or no truncation (that is the question)
% noise = randn([91,109,91]); FWHM = 3; xvals_vecs = 1:91;
% applyconvfield([50,40,60]', noise, FWHM, NaN, 0) % No truncation
% sigma = FWHM2sigma(FWHM); truncation = round(4*sigma);
% applyconvfield([50,40,60]', noise, FWHM, NaN, -1) %Default truncation
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% If a mask is supplied, mask the data
if exist('mask', 'var') && ~isnan(sum(mask(:)))
    lat_data = lat_data.*mask; %This masks the unsmoothed data
end

% Get the dimensions of the the data
Ldim = size(lat_data);

%Allow column vectors to be used as inputs
if length(Ldim) == 2 && Ldim(2) == 1
    Ldim = fliplr(Ldim); 
    lat_data = lat_data';
end

% Deal with the dimension in 1D
if Ldim(1) == 1
    Ldim = Ldim(2:end); 
end

% Get the number of dimensions
D = length(Ldim); 

%% Error checking
%--------------------------------------------------------------------------
% Error check to ensure D <= 3
if D > 3
    error('D ~= 1,2,3 has not been coded here')
end

% Ensure tval is a D by nvals matrix
if size(tval, 1) ~= D
    error(['The size of the point to evaluate must match the number of ',...
           'dimensions. I.e. the columns of the matrix must be the same',...
           'as the number of dimensions, if there is only one data point',...
           'it must be a row vector! Note that in 1D lat_data '])
end

%%  add/check optional values
%--------------------------------------------------------------------------

%% Setting the truncation and the Kernel
% Set the truncation of the kernel
if ~exist('truncation', 'var')
    truncation = -1; % setting truncation = -1 yields the default truncation below
end

% Error check to ensure truncation is the right size
if ~isequal(size(truncation),[1,1])
    error('Truncation must be a number')
end

% If Kernel is numeric take an isotropic Gaussian kernel with this FWHM
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

% Obtain the dimensions of the output
outputdim = length(Kernel(tval(:,1)));

% Initialize the matrix to store the convolution field outputs
field_vals = zeros(outputdim, size(tval, 2));

% Initialize ss 
% (used to record the sum of squares of the kernel in case you want to normalize)
ss = zeros(1, size(tval, 2));

%% Set the default xvals_vecs
% Default takes a D-dimensional cube with resolution 0 i.e. just the
% lattice point where each voxel is square with volume 1
if ~exist('xvals_vecs', 'var')
    xvals_vecs = {1:Ldim(1)}; %The other dimensions are taken case of below.
end
if ~iscell(xvals_vecs)
    xvals_vecs = {xvals_vecs};
end

% Ensure that xvals_vecs is a cell with lengh the same as the number of
% dimensions as the data
if length(xvals_vecs) < D
    increm = xvals_vecs{1}(2) - xvals_vecs{1}(1);
    for d = (length(xvals_vecs)+1):D
        xvals_vecs{d} = xvals_vecs{1}(1):increm:(xvals_vecs{1}(1)+increm*(Ldim(d)-1));
    end
end

%%  main function
%--------------------------------------------------------------------------
% If there is no truncation, obtain the indices of all the voxels on the
% lattice
if truncation == 0
    if D == 1
        xvalues_at_voxels = xvals_vecs{1}';
    elseif D == 2
        [x, y] = ndgrid(xvals_vecs{1},xvals_vecs{2});
        xvalues_at_voxels = [x(:),y(:)];
    elseif D == 3
        [x, y, z] = ndgrid(xvals_vecs{1},xvals_vecs{2},xvals_vecs{3});
        xvalues_at_voxels = [x(:),y(:),z(:)];
    end
% else
%     trunwindow = repmat(round(4*sigma), 1 ,D)'; %For non-isotropic kernel
%     we'll need to change this.
end

% Shift the data (needed if xvals_vecs is different than the default)
if truncation > 0
    % this is a fix for now, need to make it so that arbitrary xvals_vecs
    % can be used as input!! I.e. if you have irregularly spaced grid
    % points such as in MEG!
    for d = 1:D
        tval(d,:) = tval(d,:) - xvals_vecs{d}(1) + 1;
        xvals_vecs{d} = xvals_vecs{d} - xvals_vecs{d}(1) + 1;
    end
end

% Main loop
for I = 1:size(tval, 2)
    if truncation > 0
        % Obtain limits for a box around the point tval(:,I)
        lower = floor(tval(:,I) - truncation);
        upper = ceil(tval(:,I) + truncation);
        
        % For each dimension obtain the indicies (xvalues_at_voxels) of the 
        % voxels in the box specified by truncation and obtain the values
        % of the lattice data at these locations (Ytemp)
        if D == 1
            xvalues_at_voxels = (max(lower(1),  xvals_vecs{1}(1)):min(upper(1),  xvals_vecs{1}(end)))';
            Ytemp = lat_data(xvalues_at_voxels);
        elseif D == 2
            eval1 = max(lower(1),  xvals_vecs{1}(1)):min(upper(1),  xvals_vecs{1}(end));
            eval2 = max(lower(2),  xvals_vecs{2}(1)):min(upper(2),  xvals_vecs{2}(end));
            [x, y] = ndgrid(eval1,eval2);  
            xvalues_at_voxels = [x(:),y(:)];
            Ytemp = lat_data(eval1, eval2);
        elseif D == 3
            eval1 = max(lower(1),  xvals_vecs{1}(1)):min(upper(1),  xvals_vecs{1}(end));
            eval2 = max(lower(2),  xvals_vecs{2}(1)):min(upper(2),  xvals_vecs{2}(end));
            eval3 = max(lower(3),  xvals_vecs{3}(1)):min(upper(3),  xvals_vecs{3}(end));
            [x, y, z] = ndgrid(eval1,eval2,eval3);
            xvalues_at_voxels = [x(:),y(:),z(:)];
            Ytemp = lat_data(eval1, eval2, eval3);
        end
    else
        % If there is no truncation use all of the data
        Ytemp = lat_data;
    end
   
    % Obtain the kernel at all voxels in the specified truncation area
    % Need to do this for 2016 and prior compatibility
    Kernel_eval = Kernel(repmat(tval(:,I),1, length(Ytemp(:))) - xvalues_at_voxels'); 
    %     Kernel_eval = Kernel(tval(:,I) - xvalues_at_voxels'); %this is the Ith tval!
    ss(I) = sum(Kernel_eval(:).^2);
    
    % Obtain the convolution fields by multiplying the kernel evaluated at 
    % the lattice locations with the lattice data
    field_vals(:,I) = Kernel_eval*Ytemp(:);
end 

% Note could use something similar to MkRadImg to figure out which voxels are
% close to the target voxel so that you only evaluate at those voxels. i.e
% at the moment it's a square, could make it a sphere though! Ask Fabian
% what he thinks about this!

end


%             xvalues_at_voxels = (max(lower(1),  xvals_vecs{1}(1)):min(upper(1),  xvals_vecs{1}(end)))'  - xvals_vecs{1}(1) + 1;
%             eval1 = max(lower(1),  xvals_vecs{1}(1)):min(upper(1),  xvals_vecs{1}(end)) - xvals_vecs{1}(1) + 1;
%             eval2 = max(lower(2),  xvals_vecs{2}(1)):min(upper(2),  xvals_vecs{2}(end)) - xvals_vecs{2}(1) + 1;
%             eval1 = max(lower(1),  xvals_vecs{1}(1)):min(upper(1),  xvals_vecs{1}(end)) - xvals_vecs{1}(1) + 1;
%             eval2 = max(lower(2),  xvals_vecs{2}(1)):min(upper(2),  xvals_vecs{2}(end)) - xvals_vecs{2}(1) + 1;
%             eval1 = max(lower(1),  xvals_vecs{1}(1)):min(upper(1),  xvals_vecs{1}(end)) - xvals_vecs{1}(1) + 1;
%             eval2 = max(lower(2),  xvals_vecs{2}(1)):min(upper(2),  xvals_vecs{2}(end)) - xvals_vecs{2}(1) + 1;
%             eval3 = max(lower(3),  xvals_vecs{3}(1)):min(upper(3),  xvals_vecs{3}(end)) - xvals_vecs{3}(1) + 1;