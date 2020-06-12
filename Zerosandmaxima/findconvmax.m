function [maxloc, maxval] = findconvmax( lat_data, FWHM, peak_est_locs, mask, xvals_vecs, truncation)
% FINDCONVMAX( lat_data, FWHM, top_peaks ) finds the maximum in a
% convolution field by finding the largest local max.
%--------------------------------------------------------------------------
% ARGUMENTS
% lat_data      the data on a lattice to be smoothed
% FWHM          the FWHM of the kernel with which to do smoothing
% peak_est_locs     a D by npeaks matrix giving the initial estimates of the 
%               location of the peaks. If this is instead an integer: top  
%               then the top number of maxima are considered and initial 
%               locations are estimated from the underlying data. If this 
%               is not specified then it is set to 1, i.e. only considering
%               the maximum. If D = 1 and you wish to specify multiple
%               peaks rather than a number of peaks then you need to begin
%               your input as [NaN, peakestloc1, peakestloc2, ...].
%--------------------------------------------------------------------------
% OUTPUT
% % 2D
% Y = [1,1,1,1;1,2,2,1;1,2,2,1;1,1,1,1];
% [maxloc, maxval] = findconvmax(Y, 2, 1)
% fine_eval = convfield(Y, 2, 0.01, 2);
% max(fine_eval(:))
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if nargin < 4
    mask = NaN;
end
if nargin < 5
    xvals_vecs = NaN;
end
if nargin < 6
    truncation = 0;
end
[setofmaxima, setofmaxvals] = findconvpeaks(lat_data, FWHM, peak_est_locs, mask, xvals_vecs, truncation);

[~, maxind] = max(setofmaxvals);

maxloc = setofmaxima(maxind);
maxval = setofmaxvals(maxind);

end

