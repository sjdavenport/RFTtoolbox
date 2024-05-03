function LKCs = Lambda2LKC( Lambda, mask, xvalues )
% Lambda2LKCs( Lambda ) takes in a Lambda matrix corresponding to the
% variance of the derivatives of a stationary random field on the specified 
% mask and returns the LKCs.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%   Lambda:     a D by D matrix giving the variance of the derivatives
%   mask:   
%   xvalues: 
%--------------------------------------------------------------------------
% OUTPUT
%   LKCs:   a vector giving the LKCs
%--------------------------------------------------------------------------
% EXAMPLES
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Main Function Loop
%--------------------------------------------------------------------------
voxmfd  = VoxManifold( constfield( Lambda, mask, xvalues ) );
LKCs = LKC_est( voxmfd );

end

