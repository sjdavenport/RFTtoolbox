function coverage = boot_rc( data, sample_size, Kernel, resadd, mask, niters, lkc_est_version, do_spm, with_rep )
% BOOT_RC( data, sample_size, Kernel, resadd, mask, niters, lkc_est_version, do_spm, with_rep )
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  data         either an array of data or a directory containing image
%               nifti files. If a directory you need to have the BrainSTAT
%               toolbox installed in order to run this.
%  sample_size   the size of each sample to be sampled from the data
%  FWHM          the applied FWHM of the Gaussian Kernel in each direction
%         (we smooth with an istropic Kernel as is commonly done in practice)
% Optional
%  resadd       a non-negative integer giving the resolution increase.
%               Default is 1.
%  mask          a 0/1 array of size Dim which provides a mask of the data
%               the default to use is no mask i.e. 1
%  niters        the number of resamples of the data to do
%  lkc_est_version      either 'conv' or 'hpe'. Default is 'conv'
%  do_spm       additionally calculate the lkcs using SPM (i.e. under
%               stationarity
%  with_rep     0/1 whether the bootstrap samples are taken with or without
%               replacement. Default is 0 i.e. without replacement
%--------------------------------------------------------------------------
% OUTPUT
%  coverage    a structural array with entries:
%     .conv    the average coverage obtained over the number of niters
%               using the convolution approach
%     .lat     coverage obtained via evaluation on a lattice (this is
%              conservative)
%     .finelat   the coverage obtained with evaluation on the fine lattice
%              given by resadd, this is conservative but not as bad as 
%             coverage.lat. On small domains or for very smooth fields this
%             will be similar to coverage.conv, especially for high resadd. 
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'niters', 'var' )
    % default option of niters
    niters = 1000;
end

% Set the default resadd value
if ~exist('resadd', 'var')
   resadd = 1; 
end

% Set the default do_spm value
if ~exist('do_spm', 'var')
   do_spm = 1; 
end

if ~exist('lkc_est_version', 'var')
    lkc_est_version = 'conv';
elseif ~(strcmp(lkc_est_version, 'conv') || strcmp(lkc_est_version, 'hpe'))
    error('lkc_est_version must be conv or hpe');
end

% Set the default with_rep value
if ~exist('with_rep', 'var')
   with_rep = o; 
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
    spfn = @(nsubj) get_sample_fields(data, mask, nsubj, with_rep, total_nsubj, D, indexD);
elseif ischar(data)
    % Ensure that you have the BrainSTAT toolbox loaded
    if ~exist('loadsubs', 'file')
        error('You must install the BrainSTAT toolbox in order to run this option')
    end
    
    % Check the mask has the right size
    if ~isequal(size(mask), [91,109,91])
       error('The size of the mask must be 91x109x91') 
    end
    
    % Obtain the nifti files in the directory (which have either a .nii or
    % a .nii.gz extension)
    nifti_file_locs = filesindir(data, '.nii');
    
    % Obtain the total number of subjects
    total_nsubj = length(nifti_file_locs);
    
    % Determine average types (i.e. .nii or .nii.gz) of each nifti file
    use_nif = mean(nifti_type(nifti_file_locs)) - 1;
    
    % Ensure that all of the .nii/.nii.gz files are all one or the other
    if use_nif ~= 0 || use_nif ~= 1
       error('All the files in the directory must be of the same type') 
    end
    
    % Obtain the sample function
    spfn = @(nsubj) get_sample_fields_nifti(directory, nifti_file_locs, ...
                    use_nif, mask, nsubj, with_rep, total_nsubj, D, index);
else
    error('Other situations have not been coded yet')
end

% Obtain the coverage
coverage = record_coverage( spfn, sample_size, Kernel, resadd, niters, lkc_est_version, do_spm );

end

% Function to obtain a subset field
function lat_data = get_sample_fields(data, mask, nsubj, with_rep, total_nsubj, D, index)
    subset = randsample(total_nsubj,nsubj,with_rep);
    index{D+1} = subset;
    lat_data = Field(mask);
    lat_data.field = data(index{:});
end

% Function to obtain a subset field form nifti file locations
function lat_data = get_sample_fields_nifti(directory, nifti_file_locs, use_nif, mask, nsubj, with_rep, total_nsubj)
    subset = randsample(total_nsubj,nsubj,with_rep);
    lat_data = Field(mask);
    lat_data.field = loadsubs( subset, directory, use_nif, nifti_file_locs, 1 );
end
