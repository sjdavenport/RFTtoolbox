function [ curve, x ]  = ECcurve( lat_data, limits, increm, version)
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
% nvox = 100; FWHM = 5;
% lat_data = wfield(nvox);
% [ curve, x ] = ECcurve(lat_data);
% plot(x, curve)
%
% %% 2D
% Dim = [10,10]; FWHM = 5; D = length(Dim);
% lat_data = wfield(Dim);
% [ curve, x ] = ECcurve(lat_data);
% plot(x, curve)
%
% %% 2D non-Gaussian data
% Dim = [10,10]; FWHM = 5; D = length(Dim);
% lat_data = wfield(Dim, 1, 'L', 3);
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

% Obtain the curves!
for I = 1:length(x)
    x_position = sum(ECchanges <= x(I));
    if x_position == 0
        % Here before any critical points the excursion set is just 1
        % connected component.
        curve(I) = ECvals(1);
    else
        curve(I) = ECvals(x_position+1);
    end
end

end

