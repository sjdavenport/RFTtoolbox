function cut_field = cut_field( lat_data, cutoff, use_resadd )
% CUT_FIELD( field, cutoff ) takes a field and cuts off a section around
% the boundary. Useful for making stationary fields!
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  lat_data     a field whose fibersize is 1D
% Optional
%--------------------------------------------------------------------------
% OUTPUT
%  cut_field    an object of class field which is cutoff either side
%--------------------------------------------------------------------------
% EXAMPLES
% lat_data = wfield([100,1], 100); 
% FWHM = 3; resadd = 21; params = ConvFieldParams(FWHM, resadd, 0);
% smooth_field = convfield(lat_data, params);
% cf = cut_field(smooth_field, ceil(4*FWHM2sigma(FWHM)));
%
% FWHM = 3; resadd = 21; params = ConvFieldParams(FWHM, resadd, 0);
% cutoff = ceil(4*FWHM2sigma(FWHM))
% lat_data = wfield([100+2*cutoff,1], 100); 
% smooth_field = convfield(lat_data, params)
% % cutwithr = cutoff*(resadd+1);
% cf = cut_field(smooth_field, cutoff);
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

if ~exist('use_resadd', 'var')
    use_resadd = 1;
end

if use_resadd 
    if isa(lat_data, 'ConvField')
        cutoff = cutoff*(lat_data.resadd+1);
    else
        error('The use_resadd option only works if lat_data is of class ConvField');
    end
end

nminivoxels = length(lat_data.mask) - 2*cutoff;

remaining_voxels = (cutoff+1):(length(lat_data.mask)-cutoff);
cut_field = Field(lat_data.mask(remaining_voxels));
lat_spacing = lat_data.xvals{1}(2) - lat_data.xvals{1}(1);
cut_field.xvals = {(1:nminivoxels)*lat_spacing - lat_spacing + 1};
cut_field.field = lat_data.field(remaining_voxels, :);

end


% error('Not finished yet!')
% 
% % Initialize the subset field
% subset_field = Field(true(sum(lat_data.mask),1), lat_data.fibersize);
% subset_field.xvals = lat_data.xvals;
% 
% extended_mask = reshape( repmat(lat_data.mask, 1, lat_data.fibersize), [sum(lat_data.mask), lat_data.fieldsize]);

% Alternate option:
% % Set the subset field
% for I = 1:lat_data.fibersize
%     subset_field.field( = lat_data.field(extended_mask);
%     
% end

