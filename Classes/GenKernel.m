classdef GenKernel
   % GenKernel is a class to make the use of non separable kernels easy.
   %   The class implements an object consisting keeping track of the
   %   functions which need to be applied in each dimension, if a seperable
   %   kernel is used. The fields  containing 1 by D arrays. The d-th entry
   %   represents what the sperable kernel is doing in the d-th component.
   properties
      D           % the dimension of the seperable kernel
      kernel      % a 1 by D cell array containing the function handles for
                  % the kernel in each direction
      truncation  % a 1 by D array containing the value for truncation for
                  % the kernels
      dkernel     % a 1 by D cell array containing function handles for the
                  % derivatives of the kernel
      dtruncation % a 1 by D array containing the value for truncation for
                  % the dkernels
      adjust      % a 1 by D array stating whether the kernel
                  % should be shifted by adjust
   end
end