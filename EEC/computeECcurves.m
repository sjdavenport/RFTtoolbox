function [store_curves, thresholds] = computeECcurves( spfn, params, sample_size, niters )
% NEWFUN serves as a function template.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'opt1', 'var' )
   % default option of opt1
   niters = 1000;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
% Generate average EC curves
lims = [-6,6]; increm = 0.05;
thresholds = lims(1):increm:lims(2);
store_curves = zeros(niters, length(thresholds));

for J = 1:niters
    J
    lat_data = spfn(sample_size).lat_data;
    tcfield = convfield_t(lat_data, params);
    curve = ECcurve( tcfield, lims, increm);
    store_curves(J,:) = curve;
end

end

