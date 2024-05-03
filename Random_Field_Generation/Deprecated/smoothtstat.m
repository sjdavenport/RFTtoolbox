function [ tstat, setof_smoothed_fields, xbar, std_dev ] = smoothtstat( lat_data, FWHM, resadd )
% SMOOTHTSTAT( lat_data, FWHM, spacing, usespm ) takes in lattice data, and 
% returns the smooth convolution t field
%--------------------------------------------------------------------------
% ARGUMENTS
% lat_data      a Dim by nsubj array (needs nsubj > 1 in order to calculate the variance) 
% FWHM          the FWHM with which to smooth
% spacing       the interval at which to compute the convolution field.
%               I.e. if size(lat_data) = [10,20] and spacing = 0.1 the field 
%               will be evaluated at the points 1:0.1:10. Default is
%               spacing = 1.
% usespm        0/1. Whether to use spm_smooth or the inbuilt matlab
%               function convn. This is only relevant in 3D.
%--------------------------------------------------------------------------
% OUTPUT
% smoothfield   the t-statistic field that arises from smoothing the
%               individual fields and takes the t-stat.
%--------------------------------------------------------------------------
% EXAMPLES
% Dim = [5,5,5];
% nsubj = 50;
% lat_data = normrnd(0, 1, [Dim, nsubj]);
% FWHM = 1;
% smoothfield = smoothtstat(lat_data, FWHM);
% surf(smoothfield(:,:,3))
% 
% % Compare to tcfield (note quite different between spm and nonspm as
% % rough and very small data!)
% atfield = @(x) tcfield( x, lat_data, FWHM );
% smoothfield_convn = smoothtstat(lat_data, FWHM, 1, 0);
% smoothfield(3,3,3)
% smoothfield_convn(3,3,3)
% atfield([3,3,3]')
%
% % Finer comparison (perfect match up!)
% atfield = @(x) tcfield( x, lat_data, FWHM );
% smoothfield_spm = smoothtstat(lat_data, FWHM, 0.1, 1);
% smoothfield_convn = smoothtstat(lat_data, FWHM, 0.1, 0);
% point = [3,3,3]';
% latpoint = (point-1)/0.1 + 1;
% atfield(point)
% smoothfield_convn(latpoint(1), latpoint(2), latpoint(3))
% smoothfield_spm(latpoint(1), latpoint(2), latpoint(3))
%
% Dim = [5,5];
% nsubj = 50;
% lat_data = normrnd(0, 1, [Dim, nsubj])
% FWHM = 1;
% smoothfield = smoothtstat(lat_data, FWHM)
% surf(smoothfield)
%------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%------------------------------------------------------------------------
if ~exist('resadd', 'var')
    resadd = 1;
end
spacing = 1/(1+resadd);

s_lat_data = size(lat_data);
Dim = s_lat_data(1:end-1);
D = length(Dim);
Dim = (Dim-1)/spacing + 1;

setof_smoothed_fields = convfield( lat_data, FWHM, resadd, D, 0 );

if D == 1
    Dim = [1,Dim];
end
[ tstat, xbar, std_dev] = mvtstat(setof_smoothed_fields, Dim);

end

