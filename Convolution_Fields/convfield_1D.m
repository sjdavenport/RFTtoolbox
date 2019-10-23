function [cfield, xvals] = convfield_1D(data, kernel_fn, ninter, xvals_lat )
% convfield_1D(data, kernel_fn, xvals_lat, ninter )
%--------------------------------------------------------------------------
% ARGUMENTS
% data          a 1 dimensional arra  y of real values of the inital data
% kernel_fn     a function handle giving the kernel function with which to
%               convolve
% ninter        the number of subdivisions for each voxel
% max_xval      the maximum x value
%--------------------------------------------------------------------------
% OUTPUT
% cfield        a 1D array with the convolution field
% xvals         a 1D array the the x values of the voxels in the
%               convolution field
%--------------------------------------------------------------------------
% EXAMPLES
% L = 10;
% data = normrnd(0,1,1,L);
% kernel_fn = @(x) Gker(x,1);
% ninter = 100;
% cfield = convfield_1D(data, kernel_fn, ninter)
% increm = 1/ninter;
% plot(1:increm:L, cfield)
%--------------------------------------------------------------------------
% AUTHOR: Sam Davenport
if isnumeric(kernel_fn)
    kernel_fn = @(x) Gker(x,kernel_fn);
end
if nargin < 3
    ninter = 100;
end
if nargin < 4
    xvals_lat =  0:(length(data)-1);
end


if length(xvals_lat) ~= length(data)
    error('xvals_lat needs to have the same length as data');
end

if size(data,1) > 1
    error('data needs to be a row vector')
end

max_xval = max(xvals_lat);
min_xval = min(xvals_lat);
delta = (xvals_lat(2) - xvals_lat(1))/ninter; %Need to ensure the differene is always the same!

subsampled_freq = -max_xval:delta:max_xval;

data_subsampled = zeros(length(min_xval:delta:max_xval),1); %so assumes your xvalues are: 0, delta, 2delta, ..., max_freq
data_subsampled(1:ninter:length(data_subsampled)) = data;
cfield = conv(data_subsampled,kernel_fn(subsampled_freq), 'same');
xvals = min_xval:delta:max_xval;

end