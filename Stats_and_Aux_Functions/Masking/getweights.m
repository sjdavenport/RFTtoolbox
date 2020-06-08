function [ out ] = getweights( mask_hr )
% GETWEIGHTS( mask_hr )
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
% 	mask_hr     the high resolution mask
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

%%  add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'opt1', 'var' )
   % default option of opt1
   opt1 = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------

% Divide voxels in the high res mask in two in each dimension 
mask_with_divided_voxels = mask_highres( mask, 2 );

% Erode the new mask by one mask
(Use imdilate!)

% 

end

