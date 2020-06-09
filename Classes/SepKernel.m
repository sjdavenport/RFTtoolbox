classdef SepKernel
   % SepKernel is a class to make the use of seperable kernels easy.
   %   The class implements an object keeping track of the functions which
   %   need to be applied in each dimension, if a seperable kernel is used. 
   %   The fields  containing 1 by D arrays. The d-th entry represents what
   %   the sperable kernel is doing in the d-th component.
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
   methods
      function obj = SepKernel( D, FWHM )
          % SEPKERNEL( D, FWHM ) is a basic constructor for a SepKernel
          % class object. If FWHM is provided it outputs an object of class
          % Kernel representing a seperable Gaussian kernel with FWHM.
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
            obj.adjust = zeros( [ 1 D] );
          else
              if length( FWHM ) == 1 || length( FWHM ) == D
                  % create a standard SepKernel object of dimension D
                  obj = GaussKernel( SepKernel( D ), FWHM );           
              else
                error( "FWHM needs to be either of length 1 or D." )
              end
          end
      end
      
      function obj = GaussKernel( obj, FWHM )
          % GAUSSKERNEL( obj, FWHM ) constructs from a basic Kernel class
          % object a Kernel class object representing a seperable Gaussian
          % Kernel
          %----------------------------------------------------------------
          % ARGUMENTS
          % Mandatory
          %  obj   an SepKernel object
          %  FWHM  either a numeric or a numeric vector of length obj.D
          %
          %----------------------------------------------------------------
          % OUTPUT
          % obj  an object of class SepKernel representing a seperable
          %      Gaussian Kernel.
          %
          %----------------------------------------------------------------
          % EXAMPLES
          % % create an object of class SepKernel representing an isotropic
          % % Gaussian kernel
          % D = 2
          % sepK = SepKernel( D )
          % sepK = GaussKernel( sepK, 6 )
          %
          % % create an object of class SepKernel representing an seperable
          % % Gaussian kernel
          % D = 2
          % sepK = SepKernel( D )
          % sepK = GaussKernel( sepK, [ 6, 2 ] )
          %
          %----------------------------------------------------------------
          % Author: Fabian Telschow
          %----------------------------------------------------------------
          %%% check input
          % check input argument and make it a vector
          if length( FWHM ) == 1
              FWHM = FWHM * ones( [ 1 obj.D ] );
          end
          
          %%% main function
          if  length( FWHM ) == obj.D
            % fill the field kernel with Gaussians
            for d = 1:obj.D
                obj.kernel{d} = @(x) Gker( x, FWHM(d) );
            end
            
            % fill the field dkernel with derivatives of the Gaussians
            for d = 1:obj.D
                obj.dkernel{d} = @(x) Gkerderiv( x, FWHM(d) );
            end
                        
            % fill the field truncation and dtruncation
            obj.truncation  = ceil( 4 * FWHM2sigma( FWHM ) );
            obj.dtruncation = obj.truncation;
          else
              error( "FWHM needs to be either of length 1 or D." )
          end          
      end
      
      function gradobj = Gradient( obj )
          % GRADIENT( obj ) constructs a cell array containing the partial
          % derivatives of the kernel in a cell array
          %----------------------------------------------------------------
          % ARGUMENTS
          % Mandatory
          %  obj  an SepKernel object
          %
          %----------------------------------------------------------------
          % OUTPUT
          % gradobj  an SepKernel object where only the kernel
          %          specifications are included, but not the dkernel.
          %
          %----------------------------------------------------------------
          % EXAMPLES
          % % generate a separable Gaussian kernel
          % gK = SepKernel( 3, [ 2, 3, 6 ] )
          %
          % % get a cell containing the appropriate function handles in the
          % % dth entry for the dth partial derivative
          % grad_gK = Gradient( gK )
          %
          %----------------------------------------------------------------
          % Author: Fabian Telschow
          %----------------------------------------------------------------
          %%% Check mandatory input
          
          %%% main function
          % initilize the gradobj
          gradobj = SepKernel( obj.D );
          gradobj.truncation = NaN * ones( obj.D );
          
          % fill gradobj with the apropriate values
          for d = 1:obj.D
              % correct kernel functions for each dimension
              gradobj.kernel{d}    = obj.kernel;
              gradobj.kernel{d}{d} = obj.dkernel{d};
              
              % correct truncation for each dimension
              gradobj.truncation(:,d) = obj.truncation;
              gradobj.truncation(d,d) = obj.dtruncation(d);
          end
      end
      
   end
end