function grad = Gradient( obj )
% GRADIENT( obj ) constructs a cell array containing the partial
% derivatives of the kernel as Sepkernel objects.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  obj  an SepKernel object
%
%--------------------------------------------------------------------------
% OUTPUT
% grad  a 1 x obj.D cell array containing an SepKernel object for
%       each partial derivative of the inputed SepKernel.
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

%% Main function
%--------------------------------------------------------------------------

% Initialize the gradobj
grad = cell( [ 1 obj.D ] );

% Fill grad with the apropriate values
for d = 1:obj.D
  % Initialize SepKernel object for the d-th partial derivative
  % by setting it to be the obj itself
  grad{d} = obj;

  % Correct kernel functions by taking the derivative in the
  % dth component
  if isa(grad{d}.dkernel{d}, 'function_handle')
    grad{d}.kernel{d} = obj.dkernel{d};
    if iscell( obj.d2kernel )
        grad{d}.dkernel{d} = obj.d2kernel{d};
    end
  else
      [dK, d2K] = MVderiv(obj.kernel{d}, obj.D);
      grad{d}.kernel{d} = dK;
      if iscell( obj.d2kernel )
        grad{d}.dkernel{d} = d2K;
      end
  end
  
  % Correct truncation for the derivative kernel
  grad{d}.truncation(d) = obj.dtruncation(d);
end

return