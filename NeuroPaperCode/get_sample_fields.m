function spfn = get_sample_fields( data, mask )
% get_sample_fields( data, mask ) obtains a function that generates fields
% from data
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% Optional
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
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Ensure the mask is logical
mask = logical(mask);

%%  Add/check optional values
%--------------------------------------------------------------------------
% Set the default with_rep value
if ~exist('with_rep', 'var')
   with_rep = 0; 
end

%%  Main Function Loop
%--------------------------------------------------------------------------
if isnumeric(data)
    % Get mask dimensions
    [Dim, D] = getdim(mask);
    
    % Index for multidimensional coding
    indexD = repmat( {':'}, 1, D );
    
    % Obtain the size of the data
    sD = size(data);
    
    % Obtain the total number of subjects
    total_nsubj = sD(end);
    
    % Ensure the dimensions of the data and the mask match
    if ~isequal(sD(1:end-1), Dim)
        error('The dimensions of the data must match those of the mask')
    end
    
    % Obtain the sample function
    spfn = @(nsubj) get_sample_fields_mate(data, mask, nsubj, with_rep, total_nsubj, D, indexD);
 
elseif ischar(data)
    if isempty(strfind(data, '/'))
        % If the input is a folder not a file paths it defaults to the
        % featruns folder
        global featruns
        directory = [featruns, data, '_warped'];
    else
        directory = data;
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
    spfn = @(nsubj) get_sample_fields_nifti(directory, nifti_file_locs,...
                              use_nif, mask, nsubj, with_rep, total_nsubj);
else
    error('Other situations have not been coded yet')
end

end

% Function to obtain a subset field
function out = get_sample_fields_mate(data, mask, nsubj, with_rep, total_nsubj, D, index)
    out.subset = randsample(total_nsubj,nsubj,with_rep);
    index{D+1} = out.subset;
    out.lat_data = Field(mask);
    out.lat_data.field = data(index{:});
end

% Function to obtain a subset field form nifti file locations
function out = get_sample_fields_nifti(directory, nifti_file_locs, use_nif, mask, nsubj, with_rep, total_nsubj)
    out.subset = randsample(total_nsubj,nsubj,with_rep);
    out.lat_data = Field(mask); as_3D = 1;
    out.lat_data.field = loadsubs( out.subset, directory, use_nif, mask, as_3D, nifti_file_locs );
end
