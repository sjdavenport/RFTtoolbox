function spfn = getUKBspfn( RSfolder, do_gaussianize, mask )
% GETUKBSPFN( data, do_gaussianize, mask ) obtains a function that 
% generates fields from data.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  RSfolderordata    either a string giving the name of the resting state 
%                    set of data that is being worked with or a 
% Optional
%  do_gaussianize  0/1 whether or not to Gaussianize the data. Default is 1
%                  i.e. to perform Gaussianization.
%  mask         the mask for the data. If this is an array then it is taken
%               to be the mask for the data. If is it not specified then
%               the masks are obtained from the corresponding directory
%--------------------------------------------------------------------------
% OUTPUT
%  spfn         a sample function
%--------------------------------------------------------------------------
% EXAMPLES
% %% Resting State Data example
% mask = imgload('MNImask');
% spfn = getUKBspfn( 'RS_2Block', 1, mask );
% random_subset = spfn(4)
% subset = [1,3,5]; spec_subset = spfn(subset);
% subs = loadUKB([1,3,5]);
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
% Set the default do_gaussianize value
if ~exist('do_gaussianize', 'var')
    do_gaussianize = 1;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
if ischar(mask)
    % Obtain the sample function
    spfn = @(subset) get_sample_fields_nifti_intersect(subset, do_gaussianize, RSfolder);
else
    % Ensure the mask is logical
    mask = logical(mask);
    
    % Obtain the sample function
    spfn = @(subset) get_sample_fields_nifti(subset, do_gaussianize, mask, RSfolder);
end

end

% Function to obtain a subset field form nifti file locations
function lat_data = get_sample_fields_nifti(subset2use, do_gaussianize, mask, RSfolder)
    if length(subset2use) == 1
        subset2use = {subset2use};
    end
    drawset = 'cope';
    lat_data = loadUKB( subset2use, RSfolder, 'cope', drawset, mask );
    if do_gaussianize
        lat_data = Gaussianize(lat_data);
    end
end

% Function to draw a given subset with the intersection mask!
function lat_data = get_sample_fields_nifti_intersect(subset2use, do_gaussianize, RSfolder)
    if length(subset2use) == 1
        error('Not set up for random subject selection or nsubj = 1');
    end
    drawset = 'copemask';
    mask_fields = loadUKB( subset2use, RSfolder, 'mask', drawset, ones([91,109,91]));
    
    % Obtain the intersection mask by taking the product of the individual
    % subject masks
    intersection_mask = prod(mask_fields.field,4);
    
    % Obtain the lattice data (including bounding the mask
    lat_data = loadUKB( subset2use, RSfolder, 'cope', drawset, intersection_mask, RSfolder);
    if do_gaussianize
        lat_data = Gaussianize(lat_data);
    end
end

% Function to obtain a subset field form nifti file locations
% function out = get_sample_fields_nifti(subset2use, do_gaussianize)
%     if length(subset2use) == 1
%         subset2use = {subset2use};
%     end
%     [ lat_data, subj_ids, subset ] = loadUKB( subset2use, RSfolder, 'cope' );
%     if do_gaussianize
%         out.lat_data = Gaussianize(lat_data);
%     else
%         out.lat_data = lat_data;
%     end
%     out.subset = subset;
%     out.subj_ids = subj_ids;
% end

% function out = get_sample_fields_nifti_varmask(directory, nifti_file_locs, use_nif, mask_dir, nsubj, with_rep, total_nsubj)
%     out.subset = randsample(total_nsubj,nsubj,with_rep);
%     mask = ones([91,109,91]);
%     masks = loadsubs( out.subset, mask_dir, use_nif, mask, as_3D, nifti_file_locs );
%     mask = logical(prod(masks, 4));
%     [ ~, mask ] = mask_bounds( mask );
%     out.lat_data = Field(mask); as_3D = 1;
%     out.lat_data.field = loadsubs( out.subset, directory, use_nif, mask, as_3D, nifti_file_locs );
% end

% function out = get_sample_fields_nifti(directory, nifti_file_locs, use_nif, mask, nsubj, with_rep, total_nsubj)
%     out.subset = randsample(total_nsubj,nsubj,with_rep);
%     out.lat_data = Field(mask); as_3D = 1;
%     out.lat_data.field = loadsubs( out.subset, directory, use_nif, mask, as_3D, nifti_file_locs );
% end
% 
% function out = get_sample_fields_nifti_varmask(directory, nifti_file_locs, use_nif, mask_dir, nsubj, with_rep, total_nsubj)
%     out.subset = randsample(total_nsubj,nsubj,with_rep);
%     mask = ones([91,109,91]);
%     masks = loadsubs( out.subset, mask_dir, use_nif, mask, as_3D, nifti_file_locs );
%     mask = logical(prod(masks, 4));
%     [ ~, mask ] = mask_bounds( mask );
%     out.lat_data = Field(mask); as_3D = 1;
%     out.lat_data.field = loadsubs( out.subset, directory, use_nif, mask, as_3D, nifti_file_locs );
% end
