function smoothed_data = circ_conv(data, FWHM)
% circ_conv(data, FWHM) performs circular convolution of a dataset 
%--------------------------------------------------------------------------
% ARGUMENTS
% 
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% circ_conv( 1:10, 2 ) 
% circ_conv([1,1,1,0,0,0,0,0], 2)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

Kernel = @(x) Gker(x,FWHM);
kernel_param = sqrt(1/(Kernel(0)^2*2*pi)); %Why: so you don't divide by zero.
range_of_conv = round(6*kernel_param);

if length(data) > range_of_conv
    longdata = [data(end-range_of_conv+1:end), data, data(1:range_of_conv)];
else
    error('Not set up for this rubbish')
end

smoothed_data = spm_conv(longdata, FWHM);

smoothed_data = smoothed_data( (range_of_conv+1): end - range_of_conv);

end

