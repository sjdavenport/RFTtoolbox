function [smoothed_data, ss] = tor_conv(data, FWHM, D)
% tor_conv(data, FWHM) performs circular convolution of a dataset 
%--------------------------------------------------------------------------
% ARGUMENTS
% data   a Dim by nsubj vector (for some nsubj > 1)
% FWHM   the FWHM with which to smooth
%--------------------------------------------------------------------------
% OUTPUT
% smoothed_data   a vector giving the smooth data
%--------------------------------------------------------------------------
% EXAMPLES
% tor_conv( 1:10, 2 ) 
% tor_conv([1,1,1,0,0,0,0,0], 2) 
% 
% % 2D single subject
% Dim = [250,250]; white_noise = randn(Dim);
% FWHM = 25; [smoothed_data,ss] = tor_conv(white_noise, FWHM, 2);
% imagesc(smoothed_data)
%
% % 2D multiple subjects
% Dim = [250,250]; nsubj = 30; white_noise = randn([Dim, nsubj]);
% FWHM = 25; [smoothed_data,ss] = tor_conv(white_noise, FWHM, 2);
% imagesc(smoothed_data(:,:,1))
%
% % Demonstrating that the field is variance 1:
% Dim = [100,100]; nsubj = 100; white_noise = randn([Dim, nsubj]);
% FWHM = 4; [smoothed_data,ss] = tor_conv(white_noise, FWHM);
% var(smoothed_data(:)/sqrt(ss))
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

Kernel = @(x) Gker(x,FWHM);
kernel_param = sqrt(1/(Kernel(0)^2*2*pi)); %Why: so you don't divide by zero.
trunc = round(6*kernel_param);

sD = size(data);
if ~exist('D', 'var')
    D = length(sD) - 1;
else
    if length(sD) ~= D && length(sD) ~= D+1
        error('Dimension mismatch')
    end
end

% Generate the fields 
if D == 1
    extended_data = cat(1, data(end-trunc+1:end, :), data, data(1:trunc, :));
elseif D == 2
% 	other_extended_data = repmat(data,3,3); Slow way to do it!
    upperleftcorner = data(1:trunc, 1:trunc, :);
    upperrightcorner = data(1:trunc, end-trunc+1:end, :);
    lowerleftcorner = data( end-trunc+1:end, 1:trunc, :);
    lowerrightcorner = data( end-trunc+1:end, end-trunc+1:end, :);
    
    upperside = data(1:trunc, 1:end, :);
    lowerside = data(end-trunc+1:end,1:end, :);
    leftside = data(1:end, 1:trunc, :);
    rightside = data(1:end, end-trunc+1:end, :);
    
    upperrow = cat(2, lowerrightcorner, lowerside, lowerleftcorner);
    middlerow = cat(2, rightside, data, leftside);
    lowerrow = cat(2, upperrightcorner, upperside, upperleftcorner);
    
    extended_data = cat(1,upperrow,middlerow,lowerrow); 
else
    error('D >= 3 has not yet been implemented')
end

% other_smoothed_data = fconv(other_extended_data, FWHM, D);
[smoothed_data,ss] = fconv(extended_data, FWHM, D);

if D == 1
    smoothed_data = smoothed_data((trunc+1): end - trunc, :);
elseif D == 2
    smoothed_data = smoothed_data((trunc+1): end - trunc, (trunc+1): end - trunc, :);
%     other_smoothed_data = other_smoothed_data( (Dim(1) + 1):2*Dim(1), (Dim(2) + 1):2*Dim(2), :);
end

end

