function signal = square_signal( dim, radii, center_locs )
% NEWFUN
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% % EXAMPLES
% dim = [50,50]; radii = 5;
% signal = square_signal(dim, radii )
% imagesc(signal)
% signal = square_signal(dim, 4, {[25,20], [25,30]} )
% imagesc(signal)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist('center_locs', 'var')
    center_locs = {dim/2 + 1/2};
end

npeaks = length(center_locs);
if length(radii) == 1
    radii = repmat(radii, 1, npeaks);
elseif length(radii) ~= npeaks
    error('The number of peaks in radii is not the same as in centre_peaks')
end

%%  Main Function Loop
%--------------------------------------------------------------------------
signal = zeros(dim);
D = length(dim);
variable_index = cell(1,D);

for J = 1:length(center_locs)
    for I = 1:D
        variable_index{I} = floor((center_locs{J}(I)-radii(J)):(ceil(center_locs{J}(I)+radii(J))));
    end
    signal(variable_index{:}) = 1;
end

end

