function covest = rfcovest( lat_data )
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
if lat_data.D > 1
    error('Not implemented for D > 1')
end

%%  Main Function Loop
%--------------------------------------------------------------------------
covest = zeros(lat_data.masksize(1), lat_data.masksize(1));
for I = 1:lat_data.masksize(1)
    demeanedvoxdata = lat_data.field(I, :) - mean(lat_data.field(I, :));
    covest(I, :) = mean(repmat(demeanedvoxdata, lat_data.masksize(1), 1).*lat_data.field, 2);
end

end

