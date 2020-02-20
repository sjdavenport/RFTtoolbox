function [cfield, xvals, ss_vec ] = inter_conv1D( data, Kernel, increment, normalize, range_of_conv )
% inter_conv1D( data, Kernel, increment, normalize, range_of_conv ) 
% generates a convolution field.
%--------------------------------------------------------------------------
% ARGUMENTS
% data      a 1 by nvox vector 
% Kernel    a function that evaluates the kernel at each point
% increment     the constant amount that each voxel is partioned into.
% normalize     0/1 whether to normalize by the sum of squares. Default is
%               0 i.e. no normalization.
%--------------------------------------------------------------------------
% OUTPUT
% cfield        a 1D array with the convolution field
% xvals         a 1D array the the x values of the voxels in the
%               convolution field 
% ss_vec        a vector with the sum of squares at each point
%--------------------------------------------------------------------------
% EXAMPLES
% nvox = 100; xvals = 1:nvox;
% FWHM = 3;
% lat_data = normrnd(0,1,1,nvox)
% cfield = @(tval) applyconvfield(tval, lat_data, FWHM);
% plot(xvals, cfield(xvals))
% hold on
% convfield = inter_conv1D( lat_data, FWHM, 0.01);
% plot(1:0.01:nvox,convfield)
%--------------------------------------------------------------------------
% NOTES
% Potentially should change that so that (given length n) instead of 
% generating data on [1,n] it generates it on [1/2, n+1/2] which might make
% more sense in the context of the brain.
% AUTHOR: Sam Davenport.
if nargin < 3
    increment = 0.01;
end
if isnumeric(Kernel)
    if Kernel < 1
        warning('Are you sure the FWHM and increm have been written the right way around?')
    end
    Kernel = @(x) Gker(x,Kernel);
    use_gaussian = 1;
end
if nargin < 4
    normalize = 0;
end
if nargin < 5
    if use_gaussian
        kernel_param = sqrt(1/(Kernel(0)^2*2*pi)); %Why: so you don't divide by zero.
        range_of_conv = round(6*kernel_param);
    else
        error('need to set this for non-gaussian kernels')
    end
end

%Ensure that the data is a 1D row vector.
s = size(data);
if length(s) > 2
    error('Need a 1D array')
elseif s(2) < s(1)
    s = fliplr(s);
    data = data';
end
if s(1) ~= 1
    error('Need a 1D array')
end

range_of_conv = min(range_of_conv, s(2)); %ensure that the range of conv isn't greater than the image size

x = -range_of_conv:range_of_conv;

xvals = 0:increment:(length(data)-1);
cfield = zeros(1, length(xvals));

set_of_increms = 0:increment:(1-increment);
nincrems = length(set_of_increms);

ss_vec = zeros(1, length(cfield));
%The case I = 1 must be covered separately in order to include the final data point as
%well.
% kernel = exp(-(x-set_of_increms(1)).^2/(2*kernel_param^2))/sqrt(2*pi*kernel_param^2); %note set_of_increms(1) is 0 so this is symmetric so you don't need to flip it!
% kernel = Gker(x-set_of_increms(1), FWHM);
kernel_atlat = Kernel(x - set_of_increms(1)); %Even though set_of_increms(1) = 0, just including to show this is that term
smoothed_data = conv(data, kernel_atlat); %This has length length(data) + length(kernel) - 1
cfield( 1:nincrems:(nincrems*length(data))) = smoothed_data((range_of_conv + 1):end - range_of_conv); %this cuts off the first bit of incomplete overlap
ss_vec(1:nincrems:(nincrems*length(data))) = sum(kernel_atlat.^2);

for I = 2:length(set_of_increms)
%     kernel    = fliplr(exp(-(x-set_of_increms(I)).^2/(2*kernel_param^2))/(sqrt(2*pi*kernel_param^2)));
    kernel_atlat    = fliplr(Kernel(x-set_of_increms(I)));
    smoothed_data = conv(data, kernel_atlat); %This has length length(data) + length(kernel) - 1
    cfield( I:nincrems:(nincrems*(length(data)-1))) = smoothed_data((range_of_conv+1):end - range_of_conv - 1); %need to change the range here
    % have range_of_conv on one side and on the other so need to get rid of the
    % first and last lot (plus one extra as the interior elements have one
    % element less than length(data) does.
    ss_vec(I:nincrems:(nincrems*(length(data)-1))) = sum(kernel_atlat.^2);
end

if normalize == 1
    cfield(:) = cfield(:)./sqrt(ss_vec(:));
end

end
