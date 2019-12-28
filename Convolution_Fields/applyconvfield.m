function [field_vals, ss] = applyconvfield(tval, Y, Kernel, truncation, xvals_vecs)
% APPLYCONVFIELD(tval, Y, Kernel, xvals_vecs, truncation)
% calculates the convolution field at a point.
%--------------------------------------------------------------------------
% ARGUMENTS
% tval      the t values (an ndim=D by nvalues matrix) at which to evaluate 
%           the field.
% Y         the lattice field an array of size ndim.
% Kernel    the smoothing Kernel (as a function)  
% truncation    a window around the points at which to evaluate the kernel
%               setting this parameter allows for quicker computation of
%           the convolution field and has very litle effect on the values
%           of the field for kernels that have light tails such as the
%           Gaussian kernel.
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
%--------------------------------------------------------------------------
% OUTPUT
% field_vals    an nvalues length vector whose entries are the field 
%               evaluated at the locations specified by xvals
%--------------------------------------------------------------------------
% EXAMPLES  
% 1D:
% Y = [1,2,3,4];
% tval = 2; FWHM = 3;
% applyconvfield(tval, Y, FWHM)
%
% %1D Comparison with SPM smoothing:
% nvox = 100;
% xvals_vecs = 1:nvox;
% nsubj = 50;
% lat_field = normrnd(0,1,nsubj,nvox);
% field_at_voxels = applyconvfield(1:nvox, mean(lat_data'), FWHM, truncation, xvals_vecs);
% [~, ss] = applyconvfield(evalloc, lat_data(:,1)', FWHM, truncation, xvals_vecs);
% plot(field_at_voxels/sqrt(ss))
% hold on
% [smoothfield, ss_spm] = spm_conv_mod(mean(lat_data,2), FWHM);
% plot(smoothfield'/sqrt(ss_spm));
% 
% 2D:
% Y = [1,2;3,4];
% tval = [1.5,2; 3.4,2]; FWHM = 3;
% applyconvfield_gen(tval, Y, @(x) GkerMV( x, FWHM))
% applyconvfield(tval, Y, FWHM)
%
% % Different sigmas
% FWHM = 3;
% round(4*FWHM2sigma(FWHM));
%
% 3D:
% FWHM = 3;
% noise = reshape(noisegen([91,109,91], 1, 6, 1), [91,109,91]);
% tic; applyconvfield_gen([50,40,60]', noise, FWHM)
% toc
% tic; applyconvfield([50,40,60]', noise, FWHM)
% toc
%
% Truncation or no truncation (that is the question)
% noise = noisegen([91,109,91], 1, 6); FWHM = 3;
% applyconvfield([50,40,60]', noise, FWHM) % No truncation
% sigma = FWHM2sigma(FWHM); truncation = round(4*sigma);
% applyconvfield([50,40,60]', noise, FWHM, truncation)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
Ydim = size(Y);
if Ydim(1) == 1
    Ydim = Ydim(2:end);
end
D = length(Ydim);

if nargin < 4
    truncation = 0;
%     scale_var = 0;
end

if ~isequal(size(truncation),[1,1])
    error('Truncation must be a number')
end

if isnumeric(Kernel)
%     sigma = FWHM2sigma(Kernel);
    if truncation == -1
        sigma = FWHM2sigma(Kernel);
        truncation = round(4*sigma);
    end
    Kernel = @(x) GkerMV(x,Kernel);
end

if nargin < 5
    xvals_vecs = {1:Ydim(1)}; %The other dimensions are taken case of below.
end

if ~iscell(xvals_vecs)
    xvals_vecs = {xvals_vecs};
end

if length(xvals_vecs) < D
    increm = xvals_vecs{1}(2) - xvals_vecs{1}(1);
    for d = (length(xvals_vecs)+1):D
        xvals_vecs{d} = xvals_vecs{1}(1):increm:(xvals_vecs{1}(1)+increm*(Ydim(d)-1));
    end
end

if truncation == 0
    if D == 1
        xvalues_at_voxels = xvals_vecs{1}';
    elseif D == 2
        [x, y] = ndgrid(xvals_vecs{1},xvals_vecs{2});
        %     [x, y] = ndgrid(1:Ydim(1),1:Ydim(2));
        xvalues_at_voxels = [x(:),y(:)];
    elseif D == 3
        [x, y, z] = ndgrid(xvals_vecs{1},xvals_vecs{2},xvals_vecs{3});
        xvalues_at_voxels = [x(:),y(:),z(:)];
    end
% else
%     trunwindow = repmat(round(4*sigma), 1 ,D)'; %For non-isotropic kernel you'll need to change this!
end
if D > 3
    error('D ~= 1,2,3 has not been coded here') %Discuss with Armin/Tom a good way of coding this in general?
end

if size(tval, 1) ~= D
    error('The size of the point to evaluate must match the number of dimensions. I.e. the columns of the matrix must be the same as the number of dimensions, if there is only one data point it mustbe a row vector!')
end
 
outputdim = length(Kernel(tval(:,1)));
field_vals = zeros(outputdim, size(tval, 2));
ss = zeros(1, size(tval, 2));

for I = 1:size(tval, 2)
    if truncation > 0
        lower = floor(tval(:,I) - truncation);
        upper = ceil(tval(:,I) + truncation);
        if D == 1
            xvalues_at_voxels = (max(lower(1),  xvals_vecs{1}(1)):min(upper(1),  xvals_vecs{1}(end)))';
            Ytemp = Y(xvalues_at_voxels);
        elseif D == 2
            eval1 = max(lower(1),  xvals_vecs{1}(1)):min(upper(1),  xvals_vecs{1}(end));
            eval2 = max(lower(2),  xvals_vecs{2}(1)):min(upper(2),  xvals_vecs{2}(end));
            [x, y] = ndgrid(eval1,eval2);  
            xvalues_at_voxels = [x(:),y(:)];
            Ytemp = Y(eval1, eval2);
        elseif D == 3
            eval1 = max(lower(1),  xvals_vecs{1}(1)):min(upper(1),  xvals_vecs{1}(end));
            eval2 = max(lower(2),  xvals_vecs{2}(1)):min(upper(2),  xvals_vecs{2}(end));
            eval3 = max(lower(3),  xvals_vecs{3}(1)):min(upper(3),  xvals_vecs{3}(end));
            [x, y, z] = ndgrid(eval1,eval2,eval3);
            xvalues_at_voxels = [x(:),y(:),z(:)];
            Ytemp = Y(eval1, eval2, eval3);
        end
    else
        Ytemp = Y;
    end
    Kernel_eval = Kernel(repmat(tval(:,I),1, length(Ytemp(:))) - xvalues_at_voxels'); %this is the Ith tval! %Need to do this for 2016 and prior compatibility.    
%     Kernel_eval = Kernel(tval(:,I) - xvalues_at_voxels'); %this is the Ith tval!
    ss(I) = sum(Kernel_eval(:).^2);
    field_vals(:,I) = Kernel_eval*Ytemp(:);
end 

% Note can use something similar to MkRadImg to figure out which voxels are
% close to the target voxel so that you only evaluate at those voxels. i.e
% at the moment it's a square, could make it a sphere though!

end