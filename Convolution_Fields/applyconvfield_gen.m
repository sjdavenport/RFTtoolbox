function [field_vals, ss] = applyconvfield_gen(tval, Y, Kernel, xvals_vecs, scale_var)
% APPLYCONVFIELD_GEN(tval, Y, Kernel, xvals_vecs)
% calculates the locations of peaks in a random field using Newton Raphson.
%--------------------------------------------------------------------------
% ARGUMENTS
% tval      the t values (an ndim=D by nvalues matrix) at which to evaluate the field
% Y         the lattice field an array of size ndim.
% Kernel    the smoothing Kernel (as a function)  
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
% applyconvfield_gen(tval, Y, @(x) GkerMV( x, FWHM))
% 
% 2D:
% Y = [1,2;3,4];
% tval = [1.5,2; 3.4,2]; FWHM = 3;
% applyconvfield_gen(tval, Y, @(x) GkerMV( x, FWHM))
% applyconvfield_gen(tval, Y, @(x) GkerMVderiv( x, FWHM))
% applyconvfield_gen(tval, Y, @(x) GkerMVderiv2( x, FWHM))
% 
% 3D:
% FWHM = 3;
% noise = reshape(noisegen([91,109,91], 1, 6, 1), [91,109,91])
% applyconvfield_gen([50,40,60], noise, @(x) GkerMV( x, FWHM))
%--------------------------------------------------------------------------
% AUTHOR: Samuel J. Davenport
Ydim = size(Y);
if Ydim(1) == 1
    Ydim = Ydim(2:end);
    Y = reshape(Y, Ydim);
end
D = length(Ydim);

if isnumeric(Kernel)
    Kernel = @(x) GkerMV(x,Kernel);
end

if nargin < 5
    scale_var = 0;
end
if nargin < 4
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

if D == 1
    xvalues_at_voxels = xvals_vecs{1};
elseif D == 2
    [x, y] = ndgrid(xvals_vecs{1},xvals_vecs{2});
%     [x, y] = ndgrid(1:Ydim(1),1:Ydim(2));
    xvalues_at_voxels = [x(:),y(:)];
elseif D == 3
    [x, y, z] = ndgrid(xvals_vecs{1},xvals_vecs{2},xvals_vecs{3});
    xvalues_at_voxels = [x(:),y(:),z(:)];
else
    error('D ~= 1,2,3 has not been coded here') %Discuss with Armin/Tom a good way of coding this in general?
end

if size(tval, 1) ~= D
    error('The size of the point to evaluate must match the number of dimensions. I.e. the columns of the matrix must be the same as the number of dimensions, if there is only one data point it mustbe a row vector!')
end

outputdim = length(Kernel(tval(:,1)));
field_vals = zeros(outputdim, size(tval, 2));
ss = zeros(1, size(tval, 2));
for I = 1:size(tval, 2)
%     field_vals(I) = sum(Y(:).*Kernel(tval(I, :) - xvalues_at_voxels));
%     field_vals(I) = sum(Y(:)'*Kernel(tval(I, :) - xvalues_at_voxels));

%Matrix version of this:
    Kernel_eval = Kernel(repmat(tval(:,I),D) - xvalues_at_voxels'); %this is the Ith tval! %Need to do this for 2016 and prior compatibility.
%     Kernel_eval = Kernel(tval(:,I) - xvalues_at_voxels'); %this is the Ith tval!
    field_vals(:,I) = Kernel_eval*Y(:);
    ss(I) = sum(Kernel_eval(:).^2);
    if scale_var
        field_vals(:,I) = field_vals(:,I)/sqrt(ss(I));
    end
end

% Note can use something similar to MkRadImg to figure out which voxels are
% close to the target voxel so that you only evaluate at those voxels.

end