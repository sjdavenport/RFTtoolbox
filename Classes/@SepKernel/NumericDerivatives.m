function obj = NumericDerivatives( obj, h )
% NUMERICDERIVATIVES( obj ) takes Sepkernel objects and fills the dkernel
% and/or d2kernel field by numeric derivatives.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%      obj      a SepKernel object
%  Optional
%      h        stepsize for differentiation. Default is h = 0.00001
%--------------------------------------------------------------------------
% OUTPUT
% obj  an SepKernel object where the dkernel and/or d2kernel field is
%      replaced by the numeric derivative 
%
%--------------------------------------------------------------------------
% EXAMPLES
% % generate a separable Gaussian kernel
% gK = SepKernel( 3, [ 2, 3, 6 ] )
%
% % get a cell containing the appropriate function handles in the
% % dth entry for the kernel of the dth partial derivative
% grad_gK = Gradient( gK )
%
%--------------------------------------------------------------------------
% Author: Fabian Telschow
%--------------------------------------------------------------------------

%% Check mandatory input
%--------------------------------------------------------------------------
if ~exist('h', 'var')
    h = 0.00001;
end
if h < 0.00001
    warning('Small h is bad due to numerical imprecision')
end

%% Main function
%--------------------------------------------------------------------------
% Fill grad with the apropriate values
for d = 1:obj.D
  % Compute the numerical
  [dK, d2K] = MVderiv(obj.kernel{d}, 1, h);
  
  % Fill the first derivative by numeric derivative
  if isnan(obj.dkernel{d}(0)) 
    obj.dkernel{d} = dK;
  end
  
  % Fill the second derivative by numeric derivative
  if isnan(obj.d2kernel{d}(0)) 
    obj.d2kernel{d} = d2K;
  end
end

return