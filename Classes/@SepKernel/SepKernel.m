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
      function obj = SepKernel( varargin )
          % SEPKERNEL( D, FWHM ) is a basic constructor for a SepKernel
          % class object. If FWHM is provided it outputs an object of class
          % SepKernel representing a seperable Gaussian kernel with the
          % provided FWHM.
          %----------------------------------------------------------------
          % ARGUMENTS
          % Mandatory
          %  D     an numeric giving the dimension of the seperable kernel
          % Optional
          %  arg2  either a numeric or a numeric vector of length obj.D
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
          if nargin >= 1
            if isnumeric( varargin{1} )
                D =  varargin{1};
                obj.D = D;
            else
                error("First Input must be an integer giving the dimension of the underlying domain.");
            end  
          else
              error("At least the dimension of the underlying domain needs to be specified!")
          end
                      
          if nargin >= 1
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
          end
          
          if nargin >= 2
              arg2 = varargin{2};
              
              if isnumeric(arg2)
              % If second argument is numeric or 1 x D make a Gausskernel
                  if length( arg2 ) == 1 || length( arg2 ) == D
                      % create a standard SepKernel object of dimension D
                      obj = GaussKernel(D, arg2);           
                  else
                    error( "FWHM needs to be either of length 1 or D." )
                  end
              elseif isa(arg2, 'function_handle')
              % If second argument is function handle add it as kernel
              % field in all dimensions
                    for d = 1:D
                        obj.kernel{d} = arg2;
                    end
              elseif iscell(arg2)
              % If second argument is cell with function handles add it
              % in the kernel field
                for d = 1:D
                    if isa(arg2{d}, 'function_handle')
                        obj.kernel{d} = arg2{d};
                    else
                        error("The cell entries must be function handles!")
                    end
                end
              end
          end
          
          if nargin >= 3
              arg3 = varargin{3};
                            
              if isnumeric(arg3)
                  % Add the truncation field
                  if length(arg3) == 1    
                      obj.truncation   = arg3 * ones([ 1 D ]);
                      obj.dtruncation  = arg3 * ones([ 1 D ]);
                      obj.d2truncation = arg3 * ones([ 1 D ]);
                  elseif length(arg3(:)) == D
                      obj.truncation   = arg3;
                      obj.dtruncation  = arg3;
                      obj.d2truncation = arg3;
                  elseif all(size(arg3) == [3 D])
                      obj.truncation   = arg3(1,:);
                      obj.dtruncation  = arg3(2,:);
                      obj.d2truncation = arg3(3,:);                                  
                  else
                    error( "truncation input needs to be either of length 1, D or 3 x D." )
                  end
                  
              elseif isa(arg3, 'function_handle')
                    for d = 1:D
                        obj.dkernel{d} = arg3;
                    end
              elseif iscell(arg3)
                    for d = 1:D
                        if isa(arg3{d}, 'function_handle')
                            obj.dkernel{d} = arg3{d};
                        else
                            error("The cell entries must be function handles!")
                        end
                    end
              end
          end
          
          if nargin >= 4
              arg4 = varargin{4};
              
              if isnumeric(arg4)
                  % Add the truncation field
                  if length(arg4) == 1    
                      obj.adjust = arg4 * ones([ 1 D ]);
                  elseif length(arg4(:)) == D
                      obj.adjust = arg4;                      
                  else
                    error( "adjust input needs to be either of length 1, D." )
                  end          
              elseif isa(arg4, 'function_handle')
                    for d = 1:D
                        obj.d2kernel{d} = arg4;
                    end
              elseif iscell(arg4)
                    for d = 1:D
                        if isa(arg4{d}, 'function_handle')
                            obj.d2kernel{d} = arg4{d};
                        else
                            error("The cell entries must be function handles!")
                        end
                    end
              end
          end
          
          if nargin >= 5
              arg5 = varargin{5};
                            
              if isnumeric(arg5)
                  % Add the truncation field
                  if length(arg5) == 1    
                      obj.truncation   = arg5 * ones([ 1 D ]);
                      obj.dtruncation  = arg5 * ones([ 1 D ]);
                      obj.d2truncation = arg5 * ones([ 1 D ]);
                  elseif length(arg5(:)) == D
                      obj.truncation   = arg5;
                      obj.dtruncation  = arg5;
                      obj.d2truncation = arg5;
                  elseif all(size(arg5) == [3 D])
                      obj.truncation   = arg5(1,:);
                      obj.dtruncation  = arg5(2,:);
                      obj.d2truncation = arg5(3,:);                                  
                  else
                    error( "truncation input needs to be either of length 1, D or 3 x D." )
                  end
              else
                   error( "truncation input needs to be either of length 1, D or 3 x D." )
              end
          end
          
          if nargin >= 6
              arg6 = varargin{6};
              
              if isnumeric(arg6)
                  % Add the truncation field
                  if length(arg6) == 1    
                      obj.adjust = arg6 * ones([ 1 D ]);
                  elseif length(arg6(:)) == D
                      obj.adjust = arg6;                      
                  else
                    error( "adjust input needs to be either of length 1, D." )
                  end
              else
                  error( "adjust input needs to be either of length 1, D." )
              end
          end
      end
      
      % Constructs a SepKernel object from 1D Gaussian kernels
      %obj = GaussKernel( D, FWHM, adjust )
            
      % Fill derivatives by Numeric derivatives
      obj = NumericDerivatives( obj, h )
      
      % Computes the gradient for a SepKernel object
      grad = Gradient( obj )
      
      % Computes the hessian for a SepKernel object
      hess = Hessian( obj )
      
      % Evaluates a SepKernel object or its derivatives on a given grid
      vals = eval( obj, xvals, derivtype )
           
   end
end