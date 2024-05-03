function rc = store_coverage( mask, FWHM_vec, nsubj_vec, resadd, niters, filename, dirsave, use_spm )
% store_coverage( Dim, mask, FWHM_vec, nsubj_vec, use_spm, resadd, niters,dirsave )
% records and saves the coverage
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  Dim        the desired dimension of the data. If D = 1, take
%             Dim = number of voxels
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
% %% 1D coverage
% FWHM_vec = 1:0.5:6;
% nsubj_vec = 10:10:100;
% 
% Dim = 100;
% mask = true([100,1]);
% 
% store_coverage( mask, FWHM_vec, nsubj_vec, 1)
%
% %% Test
% FWHM_vec = 2;
% nsubj_vec = 10;
% 
% Dim = 100;
% mask = true([100,1]);
% store_coverage(  mask, FWHM_vec, nsubj_vec, 1, 1, 'testingjala')
% 
% %% 3D coverage
% FWHM_vec = 3:6;
% nsubj_vec = 25;
%
% Dim = [5,5,5];
% mask = true(Dim);
%
% resadd = 3; niters = 1000;
%
% store_coverage( Dim, mask, FWHM_vec, nsubj_vec, 0, resadd)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Obtain the dimensions of the data from the mask
[ rc.Dim, rc.D ] = getdim( mask );

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
    dirsave = [RFTboxloc, 'Convolution_Fields/Coverage/tstat_results/'];
end

if ~exist( 'use_spm', 'var')
    use_spm = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
rc.FWHM_vec = FWHM_vec;
rc.nsubj_vec = nsubj_vec;

if rc.D == 1
    rc.spfn = @(nsubj) temp_1D_gen(mask, rc.Dim, nsubj);
else
    rc.spfn = @(nsubj) wnfield(mask,nsubj);
end
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

if ~exist( 'filename', 'var')
    filename = [num2str(rc.D), 'D_niters_', num2str(niters), '_resadd_',num2str(resadd),'.mat'];
end
filesave = [dirsave, filename];

for I = 1:lFWHM_vec
    I
    FWHM = rc.FWHM_vec(I);
    for J = 1:lnsubj_vec
        J
        sample_size = rc.nsubj_vec(J);
        if use_spm
            coverage = record_coverage( rc.spfn, sample_size, FWHM, rc.resadd, rc.niters );
        else
            coverage = record_coverage( rc.spfn, sample_size, FWHM, rc.resadd, rc.niters, 'conv', 0 );
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

function out = temp_1D_gen(mask, Dim, nsubj)
out = wnfield(Dim, nsubj);
out.mask = mask;
end
