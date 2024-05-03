function smoothed_data = statfield( mask, FWHM, bound_mask )
% NEWFUN
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% MNImask = imgload('MNImask')
% smoothed_data = statfield( MNImask, 4 )
%--------------------------------------------------------------------------
% Copyright (C) - 2023 - Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'bound_mask', 'var' )
   % Default value
   bound_mask = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
if bound_mask == 1
    [~, mask] = mask_bounds(mask);
end
mask_padded = dilate_mask(mask, ceil(4*FWHM2sigma(FWHM))) > 0;
restriction_bounds = mask_bounds(mask_padded);

lat_data = normrnd(0,1,size(mask_padded));

[ smoothed_data, ss ] = fconv( lat_data, FWHM, length(size(mask)), 8*FWHM2sigma(FWHM));

smoothed_data = smoothed_data(restriction_bounds{:}).*mask;
smoothed_data = smoothed_data./sqrt(ss);

end

