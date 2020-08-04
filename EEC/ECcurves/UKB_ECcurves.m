function UKB_ECcurves( FWHM, nsubj_vec, RSfolder, niters )
% UKB_ECCURVES( FWHM, nsubj_vec, niters ) generates EC curves from the
% UKbiobank data
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  FWHM        the FWHM with which to smooth the datas
%  nsubj_vec   a vector giving the number of subjects you wish to use
% Optional
%  niters    the number of EC curves to calculate. Default is 1000.
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
if ~exist( 'niters', 'var' )
   % default option of opt1
   niters = 1000;
end

if ~exist('RSfolder', 'var')
    RSfolder = 'RS_2Block';
end

%%  Main Function Loop
%--------------------------------------------------------------------------
% Obtain the directory where the resting state data is saved
global featruns
directory = [featruns,RSfolder,'_warped/'];

% Obtain the location of the nifti files
nifti_file_locs = filesindir(directory, '.nii');

% Obtain the total number of subjects
total_nsubj = length(nifti_file_locs);

% Determine average types (i.e. .nii or .nii.gz) of each nifti file
use_nif = mean(nifti_type(nifti_file_locs)) - 1;

% No replacement and load in as 3D files
with_rep = 0; as_3D = 1;

% Obtain the bounded mask
mask = imgload('MNImask');
[~,mask] = mask_bounds(mask);

% Use resadd = 1 for memory
resadd = 1;

% The thresholds at which to calculate the EC curves
thresholds = -8:0.01:8;

store_curves = zeros(length(nsubj_vec),niters, length(thresholds));
for I = 1:length(nsubj_vec)
    nsubj = nsubj_vec(I);
    for J = 1:niters
        subset = randsample(total_nsubj,nsubj,with_rep);
        lat_data = Field(mask);
        lat_data.field = loadsubs( subset, directory, use_nif, mask, as_3D, nifti_file_locs );
        
        tcfield = convfield_t_Field(lat_data, FWHM, resadd);
        curve = ECcurve( tcfield, [-8,8], 0.01 );
        store_curves(I,J,:) = curve;
        save([RFTboxloc,'EEC/ECcurves/UKB_ECcurves_FWHM_', num2str(FWHM)], 'store_curves')
    end
end

end

