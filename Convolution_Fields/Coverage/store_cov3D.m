function store_cov3D( FWHM_vec, nsubj_vec, filename, niters  )
% STORE_COV3D( FWHM_vec, nsubj_vec, filename, niters  )
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  FWHM_vec   a vector specifiying the FWHM to iter over
%  nsubj_vec  a vector with numbers of subjects to iter over
% Optional
%--------------------------------------------------------------------------
% OUTPUT
% 
%--------------------------------------------------------------------------
% EXAMPLES
% store_cov3D( 1, 25, 'test', 1 )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'niters', 'var' )
   % default option of niters
   niters = 1000;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
Dim = [91,109,91];
mask = logical(imgload('MNImask'));
resadd = 1; 

tic
store_coverage( Dim, mask, FWHM_vec, nsubj_vec, resadd, niters, filename)
toc

end

