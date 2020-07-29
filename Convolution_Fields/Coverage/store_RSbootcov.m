function rc = store_RSbootcov( FWHM_vec, nsubj_vec, nifti_file_dir, filename, resadd, niters, dirsave, use_spm )
% store_coverage( Dim, mask, FWHM_vec, nsubj_vec, use_spm, resadd, niters,dirsave )
% records and saves the coverage
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  Dim        the desired dimension of the data. If D = 1, take
%             Dim = number ofvoxels
%  mask       the mask to use. If D = 1, need size(mask) = [Dim,1]
%  FWHM_vec   a vector specifiying the FWHM to iter over
%  nsubj_vec  a vector with numbers of subjects to iter over
% Optional
%  use_spm    0/1 whether to record the results of using SPM
%  resadd
%  niters
%  dirsave    the directory in which to save the coverage. If not
%             specified it uses the current directory
%--------------------------------------------------------------------------
% OUTPUT
%
%--------------------------------------------------------------------------
% EXAMPLES
% filename = 'testingruns';
% FWHM = 3; sample_size = 25;
% global feat_runs
% nifti_file_dir = [featruns,'RS_2Block_warped'];
% resadd = 1; niters = 1;
% store_RSbootcov( FWHM, sample_size, nifti_file_dir, filename, 1, niters )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'resadd', 'var' )
    % default option of resadd
    resadd = 1;
end

if ~exist( 'niters', 'var' )
    % default option of niters
    niters = 1000;
end

if ~exist( 'dirsave', 'var' )
    % default option of dirsave
    global RFTboxloc
    if isempty(RFTboxloc)
        error('You need to supply dirsave')
    else
        dirsave = [RFTboxloc, 'Convolution_Fields/Coverage/RS_results/'];
    end
end

if ~exist( 'use_spm', 'var')
    use_spm = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
mask = imgload('MNImask'); %load MNI mask (requires the BrainSTAT toolbox)

rc.FWHM_vec = FWHM_vec;
rc.nsubj_vec = nsubj_vec;

rc.resadd = resadd;
rc.niters = niters;

% Set up storage
lFWHM_vec = length(rc.FWHM_vec); lnsubj_vec = length(rc.nsubj_vec);
rc.conv = zeros(lFWHM_vec,lnsubj_vec);
rc.lat = zeros(lFWHM_vec,lnsubj_vec);
rc.finelat = zeros(lFWHM_vec,lnsubj_vec);
if use_spm
    rc.convspm = zeros(lFWHM_vec,lnsubj_vec);
    rc.latspm = zeros(lFWHM_vec,lnsubj_vec);
    rc.finelatspm = zeros(lFWHM_vec,lnsubj_vec);
end

filesave = [dirsave, filename];

for I = 1:lFWHM_vec
    I
    FWHM = rc.FWHM_vec(I);
    for J = 1:lnsubj_vec
        J
        sample_size = rc.nsubj_vec(J);
        if use_spm
            error('not implemented yet need to use boot_rc')
%             coverage = record_coverage( rc.spfn, sample_size, FWHM, rc.resadd, rc.niters, 'conv', 1 );
        else
            coverage = boot_rc( nifti_file_dir, sample_size, FWHM, mask, rc.resadd, rc.niters );
        end
        rc.conv(I,J) = coverage.conv;
        rc.lat(I,J) =  coverage.lat;
        rc.finelat(I,J) =  coverage.finelat;
        if use_spm
            rc.convspm(I,J) = coverage.convspm;
            rc.latspm(I,J) =  coverage.latspm;
            rc.finelatspm(I,J) =  coverage.finelatspm;
        end
        save(filesave, 'rc');
    end
end

end