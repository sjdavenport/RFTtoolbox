function [ out ] = convcrits( lat_data, FWHM, h )
% NEWFUN serves as a function template.
%--------------------------------------------------------------------------
% ARGUMENTS
% 
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if nargin < 3
    h = 0.00001;
end

stdfield = @(x) vf(x);

s_lat_data = size(lat_data);
D = length(s_lat_data) - 1;
nsubj = s_lat_data(end);

[ ~, setof_smoothed_fields, ~, std_dev ] = smoothtstat( lat_data, FWHM );

for I = 1:nsubj
    field = @(x) applyconvfield(x, lat_data, FWHM);
    varonefield = @(x) field(x)./stdfield(x);
    derivfield = getderivs(varonefield, D, h);

    ECn = EulerCharCrit( f, D, mask, version )
    critpoints = fzero(derivfield, initial_zero_locs);
end

function sigma = getconvvar(x)
    [~, ~, sigma] = tcfield(x, lat_data, FWHM );
end
