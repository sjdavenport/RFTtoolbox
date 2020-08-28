function [forman, kiebel, conv] = stat_smoothness_sims( mask, nsubj, FWHM_vec, niters, savefilename )
% STAT_SMOOTHNESS_SIMS runs stationary smoothness estimation simulations
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  mask       the mask on which to generate the data
%  nsubj      the number of fields to run in each simulation in order to
%             estimate the smoothness
%  FWHM_vec   a vector giving the set of FWHM with which to run the
%             simulations
% Optional
%  niters     the number of iterations to run the algorithm for. Default is
%             500.
%  savefilename   the filename (including the contains directory) in which
%                 to save the simulation results 
%--------------------------------------------------------------------------
% OUTPUT
%  three structural arrays: forman, kiebel and conv each with two options:
%   fwhm_ests  a length(FWHM_vec) by niters matrix giving the fwhm
%              estimate in each setting and iteration
%   Lambda_ests   a length(FWHM_vec) by niters matrix giving the Lambda
%              estimate in each setting and iteration
%--------------------------------------------------------------------------
% EXAMPLES
% nvox = 10; mask = true([nvox, 1]); nsubj = 10; FWHM_vec = [2,3];
% stat_smoothness_sims( mask, nsubj, FWHM_vec )
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'niters', 'var' )
   % default option of opt1
   niters = 500;
end

%%  Main Function Loop
%--------------------------------------------------------------------------

% Initialize the storage matrices for the FWHM and Lambda estimates
forman.fwhm_ests = zeros(length(FWHM_vec), niters);
forman.Lambda_ests = zeros(length(FWHM_vec), niters);
kiebel.fwhm_ests = zeros(length(FWHM_vec), niters);
kiebel.Lambda_ests = zeros(length(FWHM_vec), niters);
conv.fwhm_ests = zeros(length(FWHM_vec), niters);
conv.Lambda_ests = zeros(length(FWHM_vec), niters);

% Main loop
for I = 1:length(FWHM_vec)
    I
    FWHM = FWHM_vec(I);
    for J = 1:niters
        [forman_sim, kiebel_sim, conv_sim] = compare_smoothness_ests( mask, nsubj, FWHM );
        forman.fwhm_ests(I, J) = forman_sim.fwhm;
        forman.Lambda_ests(I, J) = forman_sim.Lambda;
        kiebel.fwhm_ests(I, J) = kiebel_sim.fwhm;
        kiebel.Lambda_ests(I, J) = kiebel_sim.Lambda;
        conv.fwhm_ests(I, J) = conv_sim.fwhm;
        conv.Lambda_ests(I, J) = conv_sim.Lambda;
    end
    
    % Save the FWHM calculations in a .mat file
    if exist('savefilename', 'var')
        save(savefilename, 'forman', 'kiebel', 'conv', 'FWHM_vec', 'niters');
    end
end

end

