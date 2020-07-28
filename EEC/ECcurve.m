function [ curve, x ]  = ECcurve( lat_data, limits, increm )
% ECcurve( lat_data, limits, ninter )
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  lat_data     a field
% Optional
%  limits       the upper and lower limits at which to calculate
%  increm       the increment
%--------------------------------------------------------------------------
% OUTPUT
% curve     the Euler chacteristic curve of the data
% x         the levels at which the Euler characteristic has been
%           calculated i.e. limits(1):increm:limits(2)
%--------------------------------------------------------------------------
% EXAMPLES
% %% 1D
% nvox = 100; FWHM = 5; D = length(Dim);
% lat_data = wnfield(nvox);
% [ curve, x ] = ECcurve(lat_data);
% plot(x, curve)
%
% %% 2D
% Dim = [10,10]; FWHM = 5; D = length(Dim);
% lat_data = wnfield(Dim);
% [ curve, x ] = ECcurve(lat_data);
% plot(x, curve)
%
% %% 3D
% Dim = [100,100,100]; FWHM = 5; D = length(Dim);
% lat_data = wnfield(Dim);
% [ curve, x ] = ECcurve(lat_data);
% plot(x, curve)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'limits', 'var' )
    % default option of limits
    limits = [-2,2];
end

if ~exist( 'ninter', 'var' )
    % default option of ninter
    increm = 0.1;
end

%%  Main Function Loop
%--------------------------------------------------------------------------

ECcalc = EulerCharCrit( lat_data.field, lat_data.D, lat_data.mask );
ECchanges =  ECcalc{1}( 2:end-1, 1 )';
ECvals = ECcalc{1}( 2:end, 2 )';

x = limits(1):increm:limits(2);
curve = zeros(1, length(x));
for I = 1:length(x)
    x_position = sum(ECchanges <= x(I));
    if x_position == 0
        % Here before any critical points the excursion set is just 1
        % connected component.
        curve(I) = 1;
    else
        curve(I) = ECvals(x_position+1);
    end
end

end

