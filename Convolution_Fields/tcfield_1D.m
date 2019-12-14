function [T, mu, sigma, d] = tcfield_1D( tval, data, xvalues_at_voxels, Kernel )
% tcfield( tval, data, xvalues_at_voxels, Kernel ) calculates a 1D-t
% convolution field given a 1D lattice with xvalues given by xvalues_at_voxels
% a vector data with the values of the lattice field at the specified
% xvalues and smoothing kernel: Kernel. This function should be used to
% evaluate a given point in t-cfield. If the whole field on a fine lattice
% is desired we would recommend using inter_conv1D for speed.
%--------------------------------------------------------------------------
% ARGUMENTS
% tval      a vector values at which to evaluate the field
% data      an nsubj by nvox matrix with each row giving the values that 
%           each lattice field (i.e. pre-smoothing) takes at the lattice points.
% xvalues_at_voxels     the x-coordinate of the lattice points. Default is
%                       to take xvalues_at_voxels = 1:length(data). I.e. to
%                   have an equally spaced initial lattice with spacing 1.
% Kernel    a function handle giving the smoothing kernel. If instead of a
%           function handle a real number is given then the kernel is set
%           to be @(x) Gker(x,Kernel) i.e. the Gaussian kernel with the
%           input parameter as it's FWHM.
%--------------------------------------------------------------------------
% OUTPUT
% T         the t-convolution field at the points given by tval
% mu        the mean estimate of the convolution fields at the points given by tval
% sigma     the standard deviation of the convolution fields at the points given by tval
% d         the cohensd-convolution field at the points given by tval
%--------------------------------------------------------------------------
% EXAMPLES
% FWHM = 4; L = 100; nsubj = 40;
% data = normrnd(0,1, nsubj, L);
% tcf = @(tval) tcfield( tval, data, 1:L, FWHM );
% plot(1:L, tcf(1:L))
% --------------------------------------------------------------------------
% AUTHOR: Samuel J. Davenport
if isnumeric(Kernel)
    if Kernel < 1
        warning('Are you sure the FWHM and increm have been written the right way around?')
    end
    Kernel = @(x) Gker(x,Kernel);
end

nsubj = size(data,1);
nevals = length(tval);
Xcfields_at_tval = zeros(nsubj, nevals);

if size(data,2) ~= length(xvalues_at_voxels)
    error('The dimensions of data and xvalues_at_voxels do not match up')
end

for I = 1:nsubj
    Xcfields_at_tval(I, :) = applyconvfield_1D(tval, xvalues_at_voxels, Kernel, data(I,:));
end

[T,mu,sigma,d] = mvtstat(Xcfields_at_tval);

end

