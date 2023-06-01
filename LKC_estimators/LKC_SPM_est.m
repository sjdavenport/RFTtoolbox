function [ L, L0 ] = LKC_SPM_est( FWHM, mask )
% LKC_SPM_est( FWHM, mask ) calculates the LKCs using SPM
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% FWHM      a vector of length the number of dimensions where the dth entry
%           is the FWHM in the dth direction. If it is just a number then
%           the kernel is assumed to be isotropic.
% mask      a D dimensional logical array giving the mask of the data
%--------------------------------------------------------------------------
% OUTPUT
%  L   a 1 x voxmfd.D vector containing the LKCs L1,..., L_voxmfd.D
%  L0  an integer containing the Euler characteristic of the mask or
%       equivalently the zeroth LKC.
%--------------------------------------------------------------------------
% EXAMPLES
% mask = ones(3,3,3);
% FWHM = [3,2,1];
% LKC_SPM_est( FWHM, mask )
%
% mask = imgload('MNImask');
% [ L, L0 ] = LKC_SPM_est( FWHM, mask )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Calculate the number of dimensions
Dim = size(mask);
D = length(Dim);

% Ensure that the mask is not logical for input into SPM
mask = double(mask);

%%  Add/check optional values
%--------------------------------------------------------------------------
if length(FWHM) < D 
    if length(FWHM) > 1
        error('Number of FWHM directions must be 1 or match the number of dimensions')
    end
   FWHM = repmat(FWHM, 1, D); 
end

%%  Main Function Loop
%--------------------------------------------------------------------------
if D == 3
    spm_resel_vec = spm_resels_vol(mask,FWHM);
elseif (D == 2 && isequal( mask, ones(Dim))) || (D == 1 && isequal( mask, ones([Dim, 1])))
    warning('2D SPM LKC esimation only works for a box atm!')
    spm_resel_vec = spm_resels(FWHM,size(mask), 'B');
else
    error('This setting has not been coded')
end

[ L, L0 ] = resel2LKC( spm_resel_vec );

end

