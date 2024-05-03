function [ xbar, std_dev, mos, tstat ] = meanmos( data, smooth_var, threeD, setnanaszero)
% MEANMOS( data ) calculates the mean and mean divided by the standard deviation
% for the dataset included in data.
%--------------------------------------------------------------------------
% ARGUMENTS
% data          An nsubj by nvox matrix with the data. Can also take in an 
%               nsubj by Dim matrix with the data but not fixed the output
%               for this yet, the output is then [1, 91, 109, 91] not
%               [91,109,91].
% smoothvar     0 or 1, specifing whether to smooth the variance (1) or not (0).
% threeD        Whether to return a 3D image or not.
%--------------------------------------------------------------------------
% OUTPUT
% xbar          The mean over subjects at each voxel.
% std           The estimate of the standard deviation at each voxel.
% mos           xbar/std.
% tstat         sqrt(nsubj)*mos: the t-statistic at each voxel.
%--------------------------------------------------------------------------
% EXAMPLES
% noise = noisegen([91,109,91], 20, 6, 1);
% [xbar, std, ~, tstat] = meanmos(noise);
% size(xbar)
%
% nsubj = 20;
% noise = noisegen([91,109,91], nsubj, 6, 3);
% [xbar, std_dev, ~, tstat] = meanmos(noise);
% vox = 150;
% noise_at_vox = noise(:, vox);
% mu = mean(noise_at_vox)
% sigmatilde = std(noise_at_vox)
% sqrt(nsubj)*mu/sigmatilde
% tstat(vox)
%--------------------------------------------------------------------------
% SEE ALSO
% mat
if nargin < 2
    smooth_var = 0;
end
if nargin < 3
    threeD = 0;
end
if nargin < 4
   setnanaszero = 1;
end
sD = size(data);
nSubj = sD(1);

% if sD(end) > sD(1)
%     warning('remember for now meanmos requires that the first corrordinate is nSubj')
% end

xbar = mean(data);
sq_xbar = mean(data.^2);
    
est_var = (nSubj/(nSubj-1))*(sq_xbar - (xbar.^2));
if smooth_var == 1
    est_var = smooth3Dvect( est_var, 3 )';
end

std_dev = sqrt(est_var);
mos = xbar./std_dev;

if threeD
    global stdsize
    mos = reshape(mos, stdsize);    
    xbar = reshape(xbar, stdsize);
    std_dev = reshape(std_dev, stdsize);
end

%NEW LINE!!! May need some testing!
if setnanaszero
    mos(isnan(mos))=0;
end

tstat = mos*sqrt(nSubj);

end

