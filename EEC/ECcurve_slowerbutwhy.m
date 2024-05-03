function [ curve, x ]  = ECcurve_alternatve( lat_data, limits, increm, version)
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

if ~exist( 'increm', 'var' )
    % default option of ninter
    increm = 0.1;
end

if ~exist( 'version', 'var' )
    % default option of ninter
    version = "C";
end

%%  Main Function Loop
%--------------------------------------------------------------------------

% Calculate the points at which the obsered EC changes
ECcalc = EulerCharCrit( lat_data.field, lat_data.D, lat_data.mask, version );
ECchanges =  ECcalc{1}( 2:end-1, 1 )';

% Findt he values of the EC after those changes
ECvals = ECcalc{1}( 2:end, 2 )';

% Obtain the vector at which to evaluate the curve
x = limits(1):increm:limits(2);

% Initialize the storage curve
curve = zeros(1, length(x));

% Obtain the curves, this is slow and can be improved by setting all values
% between change areas to the correct value.
xlocofchange_new = 1;
for I = 1:length(ECchanges)
    xlocofchange_old = xlocofchange_new;
    xlocofchange_new = find(abs(x - ECchanges(I)) < increm/2, 1, 'last');
    curve(xlocofchange_old:xlocofchange_new) = ECvals(I);
end

end

