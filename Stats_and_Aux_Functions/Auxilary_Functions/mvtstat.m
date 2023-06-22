function [tstat, xbar, std_dev, cohensd] = mvtstat( data, Dim, nansaszeros )
% MVTSTAT( data, threeD, nansaszeros ) computes the multivariate t-statistic.
%--------------------------------------------------------------------------
% ARGUMENTS
% data          a Dim by nsubj array with the data.
% Dim           the dimension to return.
% nansaszeros   0/1 whether to return NaN entries as zeros. Default is 0.
%--------------------------------------------------------------------------
% OUTPUT
% tstat         the one sample t-statistic at each voxel.
% xbar
% std_dev       the standard deviation at each voxel (calculated using the
%               unbiased estimate of the variance)
%--------------------------------------------------------------------------
% EXAMPLES
% mvtstat(normrnd(0,1,1,100))
%
% noise = noisegen([91,109,91], 20, 6, 0);
% tstat = mvtstat(noise);
%
% Dim = [100,100];
% noise = noisegen(Dim, 20, 2, 0);
% tstat = mvtstat(noise, Dim);
%
% Dim = 100;
% noise = noisegen(Dim, 20, 2, 0);
% tstat = mvtstat(noise);
% 
% nsubj = 20;
% noise = noisegen([91,109,91], nsubj, 6, 3);
% tstat = mvtstat(noise);
% vox = 150;
% noise_at_vox = noise(:, vox);
% muuuu = mean(noise_at_vox)
% sigmatilde = std(noise_at_vox)
% sqrt(nsubj)*mu/sigmatilde
% tstat(vox)
%
% % For comparison to python code:
% a = reshape(1:12,4,3)
% [tstat, xbar, std_dev] = mvtstat(data)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
sD = size(data);

if nargin < 2
    Dim = sD(1:(end-1));
end
if nargin < 3
    nansaszeros = 0;
end

D = length(sD) - 1;
nsubj = sD(end);

xbar = mean(data, (D+1));
sq_xbar = mean(data.^2, (D+1));
    
est_var = (nsubj/(nsubj-1))*(sq_xbar - (xbar.^2)); %This is the population estimate!
std_dev = real(sqrt(est_var)); % Some matlab bug sometimes makes this imaginary! See below for commented example of the bug if the real is removed
% std_dev = sqrt(est_var);

if isequal(prod(Dim),  prod(sD(1:end-1))) && (length(Dim) > 1)
    xbar = reshape(xbar, Dim);
    std_dev = reshape(std_dev, Dim);
else
    xbar = xbar(:);
    std_dev = std_dev(:);
end

tstat = sqrt(nsubj)*xbar./std_dev;
cohensd = xbar./std_dev;

if nansaszeros
    tstat = nan2zero(tstat);
    cohensd = nan2zero(cohensd);
end

end

% Example of the weird imaginary value error that can occasionally occur.
% dim = [100,100];
% nsubj = 100;
% effectsize = 0.5;
% field_type = 'L';
% field_params = 3;
% lat_data = wfield(dim, nsubj, field_type, field_params);
% 
% signal = peakgen( effectsize, 8, 10, dim, {[20,20], [20,80], [50,50], [80,20], [80,80]});
% lat_data.field = lat_data.field + signal;
% 
% FWHM = 10;
% tstat_orig = convfield_t(lat_data, FWHM).field;
% pvals_orig = tstat_pval(tstat_orig, nsubj-1, 0);
% smooth_fields = convfield(lat_data, FWHM);
% smooth_fields_gauss = Gaussianize(smooth_fields.field);
% % Need to investigate why occasionally get imaginary parts!
% 
% % smooth_fields_gauss = sqrt(smooth_fields.field);
% tstat_gauss = mvtstat(smooth_fields_gauss.field);
% pvals_gauss = tstat_pval(tstat_gauss, nsubj-1, 0);
