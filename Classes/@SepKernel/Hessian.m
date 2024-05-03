function hess = Hessian( obj ) 
% Hessian( obj ) constructs a cell array containing the second
% partial derivatives of the kernel as Sepkernel objects.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  obj  an SepKernel object
%
%--------------------------------------------------------------------------
% OUTPUT
% hess  a obj.D x obj.D cell array containing an SepKernel object
%       for each element of the Hessian matrix of the inputed
%       SepKernel.
%
%--------------------------------------------------------------------------
% EXAMPLES
% % generate a separable Gaussian kernel
% gK = SepKernel( 3, [ 2, 3, 6 ] )
%
% % get a cell containing the appropriate function handles in the
% % dth entry for the dth partial derivative
% hessian_gK = Hessian( gK )
%
%--------------------------------------------------------------------------
% Author: Fabian Telschow
%--------------------------------------------------------------------------

%% Check mandatory input
%--------------------------------------------------------------------------

%% Main function
%--------------------------------------------------------------------------

% Initialize the gradobj
hess = cell( [ obj.D obj.D ] );

% Fill hess with the appropriate values
for d = 1:obj.D
  for dd = 1:obj.D
      % Initialize SepKernel object for the d-th partial
      % derivative by setting relevant parts to the obj itself
      hess{d,dd} = SepKernel( obj.D );
      hess{d,dd}.kernel = obj.kernel;
      hess{d,dd}.truncation = obj.truncation;

      if d == dd
          % Diagonals contain second derivatives
          hess{d,dd}.kernel{d} = obj.d2kernel{d};
          % Correct truncation for the derivative kernel
          hess{d,dd}.truncation(d) = obj.d2truncation(d);
      else
          % Correct kernel functions by taking the derivative
          % in the d-th component and the dd-th
          hess{d,dd}.kernel{d} = obj.dkernel{d};
          hess{d,dd}.kernel{dd} = obj.dkernel{dd};

          % Correct truncation for the derivative kernel
          hess{d,dd}.truncation(d) = obj.dtruncation(d);
          hess{d,dd}.truncation(dd) = obj.dtruncation(dd);
      end
  end
end

return