function [tstat, xbar, std_dev, cohensd] = mvtstat( data, Dim, nansaszeros )
% MVTSTAT( data, threeD, nansaszeros ) computes the multivariate t-statistic.
%--------------------------------------------------------------------------
% ARGUMENTS
% data          A D by nsubj array with the data.
% Dim           the dimension to return.
% nansaszeros   0/1 whether to return NaN entries as zeros. Default is 0.
%--------------------------------------------------------------------------
% OUTPUT
% tstat         the one sample t-statistic at each voxel.
%--------------------------------------------------------------------------
% EXAMPLES
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
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport.
% warning('Have changed nsubj by D to D by nsubj, so may cause errors!')
if nargin < 2
    Dim = NaN;
end
if nargin < 3
    nansaszeros = 0;
end

sD = size(data);
D = length(sD) - 1;
nsubj = sD(end);

xbar = mean(data, (D+1));
sq_xbar = mean(data.^2, (D+1));
    
est_var = (nsubj/(nsubj-1))*(sq_xbar - (xbar.^2)); %This is the population estimate!
std_dev = sqrt(est_var);

if Dim == 1
    stdsize = [91,109,91];
    xbar = reshape(xbar, stdsize);
    std_dev = reshape(std_dev, stdsize);
elseif isequal(prod(Dim),  prod(sD(1:end-1)))
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

