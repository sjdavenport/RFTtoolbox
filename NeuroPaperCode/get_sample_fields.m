function spfn = get_sample_fields( RSfolder, mask, do_gaussianize )
% get_sample_fields( data, mask ) obtains a function that generates fields
% from data
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  RSfolder     the name of the resting state set of data that is being
%               worked with 
% Optional
%  mask         the mask for the data. If this is an array then it is taken
%               to be the mask for the data. If is it not specified then
%               the masks are obtained from the corresponding directory
% do_gaussianize    0/1 should you Gaussianize the data or not. Default is
%                   0 i.e. no Gaussianization
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% %% %% 1D Examples
% %% Simple 1D example
% tot_nsubj = 30; nvox = 100; data = normrnd(0,1,nvox, tot_nsubj); 
% mask = true(1,nvox)';
% spfn = get_sample_fields( data, mask );
% 
% %% %% 2D Examples
% %% Simple 2D example
% tot_nsubj = 50; Dim = [20,30]; data = normrnd(0,1,[Dim,tot_nsubj]); 
% mask = true(Dim);
% spfn = get_sample_fields( data, mask );
% spfn(20).lat_data
%
% %% Resting State Data example
% mask = imgload('MNImask');
% spfn = get_sample_fields( 'RS_2Block', mask );
% a = spfn(20)
%
% directory =
% '/vols/Scratch/ukbiobank/nichols/SelectiveInf/feat_runs/RS_2Block_warped/'
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Ensure the mask is logical
mask = logical(mask);
D = length(size(mask));

%%  Add/check optional values
%--------------------------------------------------------------------------
% Set the default with_rep value
if ~exist('with_rep', 'var')
   with_rep = 0; 
end

if ~exist('do_gaussianize', 'var')
    do_gaussianize = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
if isnumeric(RSfolder)
    % Index for multidimensional coding
    indexD = repmat( {':'}, 1, D );
    
    % Obtain the size of the data
    sD = size(RSfolder);
    
    % Obtain the total number of subjects
    total_nsubj = sD(end);
    
    % Obtain the sample function
    if isequal(size(RSfolder), size(mask))
        spfn = @(nsubj) get_sample_fields_mate_varmask(RSfolder, mask, nsubj, with_rep, total_nsubj, D, indexD);
    else
        spfn = @(nsubj) get_sample_fields_mate(RSfolder, mask, nsubj, with_rep, total_nsubj, D, indexD, do_gaussianize);
    end
 
elseif ischar(RSfolder)
    if isempty(strfind(RSfolder, '/'))
        % If the input is a folder not a file paths it defaults to the
        % featruns folder
        global featruns
        directory = [featruns, RSfolder, '_warped/'];
    else
        directory = RSfolder;
    end
    
    if ~exist('mask', 'var')
        mask_dir = [featruns, RSfolder, '_masks_warped/'];
    else
        mask_dir = 0;
    end
    
    % Obtain the nifti files in the directory (which have either a .nii or
    % a .nii.gz extension)
    nifti_file_locs = filesindir(directory, '.nii');
    
    % Obtain the total number of subjects
    total_nsubj = length(nifti_file_locs);
    
    % Determine average types (i.e. .nii or .nii.gz) of each nifti file
    use_nif = mean(nifti_type(nifti_file_locs)) - 1;
    
    % Ensure that all of the .nii/.nii.gz files are all one or the other
    if (use_nif ~= 0) && (use_nif ~= 1)
       error('All the files in the directory must be of the same type') 
    end
    
    % Obtain the sample function
    if isstr(mask_dir)
        spfn = @(nsubj) get_sample_fields_nifti_varmask(directory, ...
             nifti_file_locs, use_nif, mask_dir, nsubj, with_rep, total_nsubj);
    else
        [ ~, mask ] = mask_bounds( mask );
        spfn = @(nsubj) get_sample_fields_nifti(directory, nifti_file_locs,...
                              use_nif, mask, nsubj, with_rep, total_nsubj);
    end
else
    error('Other situations have not been coded yet')
end

end
% (RSfolder, mask, nsubj, with_rep, total_nsubj, D, indexD, do_gaussianize)
% Function to obtain a subset field
function out = get_sample_fields_mate(data, mask, nsubj, with_rep, total_nsubj, D, index, do_gaussianize)
    out.subset = randsample(total_nsubj,nsubj,with_rep);
    index{D+1} = out.subset;
    out.lat_data = Field(mask);
    out.lat_data.field = data(index{:});
    if do_gaussianize && nsubj > 1
        out.lat_data = Gaussianize(out.lat_data);
    end
end

% Function to obtain a subset field with variable masks
function out = get_sample_fields_mate_varmask(data, mask_data, nsubj, with_rep, total_nsubj, D, index)
    out.subset = randsample(total_nsubj,nsubj,with_rep);
    index{D+1} = out.subset;
    mask = logical(prod(mask_data(index{:}), D+1));
    out.lat_data = Field(mask);
    out.lat_data.field = data(index{:});
end

% Function to obtain a subset field form nifti file locations
function out = get_sample_fields_nifti(directory, nifti_file_locs, use_nif, mask, nsubj, with_rep, total_nsubj)
    out.subset = randsample(total_nsubj,nsubj,with_rep);
    out.lat_data = Field(mask); as_3D = 1;
    out.lat_data.field = loadsubs( out.subset, directory, use_nif, mask, as_3D, nifti_file_locs );
end

function out = get_sample_fields_nifti_varmask(directory, nifti_file_locs, use_nif, mask_dir, nsubj, with_rep, total_nsubj)
    out.subset = randsample(total_nsubj,nsubj,with_rep);
    mask = ones([91,109,91]);
    masks = loadsubs( out.subset, mask_dir, use_nif, mask, as_3D, nifti_file_locs );
    mask = logical(prod(masks, 4));
    [ ~, mask ] = mask_bounds( mask );
    out.lat_data = Field(mask); as_3D = 1;
    out.lat_data.field = loadsubs( out.subset, directory, use_nif, mask, as_3D, nifti_file_locs );
end
