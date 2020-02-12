function [ tstat, xbar, setof_smoothed_fields ] = smoothtstat( lat_data, FWHM )
% SMOOTHTSTAT( lat_data ) takes in lattice data, smoothes it with a certain
% FWHM and returns the t-statistic evaluated on the lattice.
%--------------------------------------------------------------------------
% ARGUMENTS
% lat_data      a Dim by nsubj array 
% FWHM          the FWHM with which to smooth
%--------------------------------------------------------------------------
% OUTPUT
% smoothfield   the t-statistic field that arises from smoothing the
%               individual fields and takes the t-stat.
%--------------------------------------------------------------------------
% EXAMPLES
% Dim = [5,5,5];
% nsubj = 50;
% lat_data = normrnd(0, 1, [Dim, nsubj])
% FWHM = 1;
% smoothfield = smoothtstat(lat_data, FWHM)
% surf(smoothfield(:,:,3))
%
% Dim = [5,5];
% nsubj = 50;
% lat_data = normrnd(0, 1, [Dim, nsubj])
% FWHM = 1;
% smoothfield = smoothtstat(lat_data, FWHM)
% surf(smoothfield)
%------------------------------------------------------------------------
% AUTHOR: Samuel Davenport

s_lat_data = size(lat_data);
nsubj = s_lat_data(end);
Dim = s_lat_data(1:end-1);
D = length(Dim);
if D == 1
    Dim = [1,Dim];
end

setof_smoothed_fields = zeros(s_lat_data);

if D == 3
    smoothfield_subj = zeros(Dim);
    for I = 1:nsubj
        spm_smooth(lat_data(:,:,:,I), smoothfield_subj, FWHM);
        setof_smoothed_fields(:,:,:,I) = smoothfield_subj;
    end
elseif D == 2
    for I = 1:nsubj
        setof_smoothed_fields(:,:,I) = spm_conv(lat_data(:,:,I), FWHM);
    end
elseif D == 1
    for I = 1:nsubj
        setof_smoothed_fields(:,I) = spm_conv(lat_data(:,I), FWHM);
    end
else
    error('Not implemented for D > 3')
end

[ tstat, xbar ] = mvtstat(setof_smoothed_fields, Dim);

end

