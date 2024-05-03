function [ curve, u ]  = ECcurve( field, u, version)
% ECcurve( lat_data, x, version ) evaluates the EC curve of the given field
% at the locations of the vector x.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  field   an object of Class Field
% Optional
%  u      a vector giving the levels where the EC curve is evaluated.
%         Default: u = -3:0.1:3;
%--------------------------------------------------------------------------
% OUTPUT
% curve     an object of Class field containing the Euler chacteristic
%           curves of the input field at the levels giving by u
%--------------------------------------------------------------------------
% EXAMPLES
% u = -5:0.1:5;
% %% 1D
% nvox = 100; FWHM = 5;
% lat_data = wfield(nvox, 10);
% [ curve, u ] = ECcurve(lat_data);
% plot(curve)
% 
% %% 2D
% Dim = [10,10]; FWHM = 5; D = length(Dim);
% lat_data = wfield(Dim, 10);
% [ curve, u ] = ECcurve(lat_data);
% plot(curve)
% 
% %% 2D non-Gaussian data
% Dim = [10,10]; FWHM = 5; D = length(Dim);
% lat_data = wfield(Dim, 1, 'L', 3);
% [ curve, u ] = ECcurve(lat_data);
% plot(curve)
% 
% %% 3D
% Dim = [100,100,100]; FWHM = 5; D = length(Dim);
% lat_data = wfield(Dim, 10);
% [ curve, u ] = ECcurve(lat_data);
% plot(curve)
%--------------------------------------------------------------------------
% AUTHOR: Fabian Telschow, Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'u', 'var' )
    % default option of limits
    u = -3:0.1:3;
end

if ~exist( 'version', 'var' )
    % default option of ninter
    version = "C";
end

%%  Main Function Loop
%--------------------------------------------------------------------------

% Calculate the points at which the observed EC changes
ECcalc = EulerCharCrit( field.field, field.D, field.mask, version );

% Initialize the storage curve
curve = zeros(length(u), length(ECcalc));

for l = 1:length(ECcalc)
    % Get critical values where the ECchanges
    ECchanges =  ECcalc{l}( 2:end-1, 1 )';

    % Find the change of the EC at the critical values
    ECvals = ECcalc{l}( 2:end, 2 )';

    % Obtain the curves!
    for I = 1:length(u)
        x_position = sum(ECchanges <= u(I));
        if x_position == 0
            % Here before any critical points the excursion set is just 1
            % connected component.
            curve(I,l) = ECvals(1);
        else
            curve(I,l) = ECvals(x_position+1);
        end
    end
end

curve = Field(curve, 1, {u});
end

