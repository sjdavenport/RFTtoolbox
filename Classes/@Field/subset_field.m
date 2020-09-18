function [ out ] = subset_field( in )
% SUBSET_FIELD( lat_data ) returns a 2 dimensional field consists of the
% voxels in the mask and the data at those voxels.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  lat_data     a field whose fibersize is 1D
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

error('Not finished yet!')

% Initialize the subset field
subset_field = Field(true(sum(lat_data.mask),1), lat_data.fibersize);
subset_field.xvals = lat_data.xvals;

extended_mask = reshape( repmat(lat_data.mask, 1, lat_data.fibersize), [sum(lat_data.mask), lat_data.fieldsize]);

end

% Alternate option:
% % Set the subset field
% for I = 1:lat_data.fibersize
%     subset_field.field( = lat_data.field(extended_mask);
%     
% end

