classdef SepKernel
   % SepKernel is a class to make the use of seperable kernels easy.
   %   The class implements an object keeping track of the functions which
   %   need to be applied in each dimension, if a seperable kernel is used. 
   %   The fields  containing 1 by D arrays. The d-th entry represents what
   %   the sperable kernel is doing in the d-th component.
   properties (Access = private)
      D           {mustBePositive, mustBeInteger} % the dimension of the seperable kernel
   end
   properties
      kernel      cell % a 1 by D cell array containing the function handles for the kernel in each direction
      truncation  % a 1 by D array containing the value for truncation for the kernels
      dkernel     cell % a 1 by D cell array containing function handles for the derivatives of the kernel
      dtruncation % a 1 by D array containing the value for truncation for the dkernels
      d2kernel    cell  % a 1 by D cell array containing function handles for the second derivatives of the kernel
      d2truncation % a 1 by D array containing the value for truncation for the d2kernels
      adjust       % a 1 by D array stating whether the kernel should be shifted by adjust
   end
   methods
      function obj = SepKernel( D, FWHM )
          % SEPKERNEL( D, FWHM ) is a basic constructor for a SepKernel
          % class object. If FWHM is provided it outputs an object of class
          % SepKernel representing a seperable Gaussian kernel with the
          % provided FWHM.
          %----------------------------------------------------------------
          % ARGUMENTS
          % Mandatory
          %  D     an numeric giving the dimension of the seperable kernel
          % Optional
          %  FWHM  either a numeric or a numeric vector of length obj.D
          %
          %----------------------------------------------------------------
          % OUTPUT
          % obj  an object of class SepKernel where the fields are
          %      initialised. If FWHM is provided obj is an object of class
          %      SepKernel representing a seperable Gaussian Kernel.
          %
          %----------------------------------------------------------------
          % EXAMPLES
          % % create a standard kernel object
          % D = 2
          % sepK = SepKernel( D )
          %
          % % create a Gaussian kernel object with FWHM
          % D = 3;
          % FWHM = [3 2 6];
          % gaussK = SepKernel( D, FWHM )
          %
          % % create an isotropic Gaussian kernel object
          % D = 3;
          % gaussK = SepKernel( D, 6 )
          %
          %----------------------------------------------------------------
          % Author: Fabian Telschow
          %----------------------------------------------------------------
          if nargin == 1
            obj.D = D;
            
            obj.kernel = cell( [ 1 D ] );
            for d = 1:D
                obj.kernel{d} = @(x) NaN;
            end
            obj.truncation = NaN * ones( [ 1 D ] );
            
            obj.dkernel = cell( [ 1 D ] );
            for d = 1:D
                obj.dkernel{d} = @(x) NaN;
            end
            obj.dtruncation = NaN * ones( [ 1 D ] );

            obj.d2kernel = cell( [ 1 D ] );
            for d = 1:D
                obj.d2kernel{d} = @(x) NaN;
            end
            obj.d2truncation = NaN * ones( [ 1 D ] );
            
            obj.adjust = zeros( [ 1 D] );
            
          else
              if length( FWHM ) == 1 || length( FWHM ) == D
                  % Create a standard SepKernel object of dimension D
                  obj = GaussKernel( SepKernel( D ), FWHM );           
              else
                error( "FWHM needs to be either of length 1 or D." )
              end
          end
      end
      
      % Constructs a SepKernel object from 1D Gaussian kernels
      obj = GaussKernel( obj, FWHM );
      
      % Computes the gradient for a SepKernel object
      grad = Gradient( obj )
      
      % Computes the hessian for a SepKernel object
      hess = Hessian( obj )
      
      % Evaluates a SepKernel object or its derivatives on a given grid
      vals = eval( obj, xvals, derivtype )
           
   end
end