classdef SepKernel
   % SepKernel is a class to make the use of seperable kernels easy.
   %   The class implements an object keeping track of the functions which
   %   need to be applied in each dimension, if a seperable kernel is used. 
   %   The fields  containing 1 by D arrays. The d-th entry represents what
   %   the sperable kernel is doing in the d-th component.
   properties (Access = private)
      D           {mustBePositive} % the dimension of the seperable kernel
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
          %%% Check input
          % Check input argument and make it a vector
          if length( FWHM ) == 1
              FWHM = FWHM * ones( [ 1 obj.D ] );
          end
          
          %%% Main function
          if  length( FWHM ) == obj.D
            % Fill the field kernel with Gaussians
            for d = 1:obj.D
                obj.kernel{d} = @(x) Gker( x, FWHM(d) );
            end
            
            % Fill the field dkernel with derivatives of the Gaussians
            for d = 1:obj.D
                obj.dkernel{d} = @(x) Gkerderiv( x, FWHM(d) );
                obj.d2kernel{d} = @(x) Gkerderiv2( x, FWHM(d) );
            end
                        
            % Fill the field truncation and dtruncation
            obj.truncation   = ceil( 4 * FWHM2sigma( FWHM ) );
            obj.dtruncation  = obj.truncation;
            obj.d2truncation = obj.truncation;
          else
              error( "FWHM needs to be either of length 1 or D." )
          end
      end
      
      function grad = Gradient( obj )
          % GRADIENT( obj ) constructs a cell array containing the partial
          % derivatives of the kernel as Sepkernel objects.
          %----------------------------------------------------------------
          % ARGUMENTS
          % Mandatory
          %  obj  an SepKernel object
          %
          %----------------------------------------------------------------
          % OUTPUT
          % grad  a 1 x obj.D cell array containing an SepKernel object for
          %       each partial derivative of the inputed SepKernel.
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
          
          %%% Main function
          % Initialize the gradobj
          grad = cell( [ 1 obj.D ] );
          
          % Fill grad with the apropriate values
          for d = 1:obj.D
              % Initialize SepKernel object for the d-th partial derivative
              % by setting it to be the obj itself
              grad{d} = obj;
              
              % Correct kernel functions by taking the derivative in the
              % dth component
              grad{d}.kernel{d} = obj.dkernel{d};
              if iscell( obj.d2kernel )
                grad{d}.dkernel{d} = obj.d2kernel{d};
              end
              % Correct truncation for the derivative kernel
              grad{d}.truncation(d) = obj.dtruncation(d);
          end
      end
      
      function hess = Hessian( obj )
          % Hessian( obj ) constructs a cell array containing the second
          % partial derivatives of the kernel as Sepkernel objects.
          %----------------------------------------------------------------
          % ARGUMENTS
          % Mandatory
          %  obj  an SepKernel object
          %
          %----------------------------------------------------------------
          % OUTPUT
          % hess  a obj.D x obj.D cell array containing an SepKernel object
          %       for each element of the Hessian matrix of the inputed
          %       SepKernel.
          %
          %----------------------------------------------------------------
          % EXAMPLES
          % % generate a separable Gaussian kernel
          % gK = SepKernel( 3, [ 2, 3, 6 ] )
          %
          % % get a cell containing the appropriate function handles in the
          % % dth entry for the dth partial derivative
          % hessian_gK = Hessian( gK )
          %
          %----------------------------------------------------------------
          % Author: Fabian Telschow
          %----------------------------------------------------------------
          %%% Check mandatory input
          
          %%% Main function
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
      end
      
      function vals = eval( obj, xvals, derivtype )
          % eval( obj ) evaluates a SepKernel object or its derivatives on
          % a given grid.
          %----------------------------------------------------------------
          % ARGUMENTS
          % Mandatory
          %  obj   an SepKernel object
          %  xvals an 1 by D cell array containing on which grid values to
          %        evaluate the d-th kernel
          % Optional
          %  derivtype 0/1/2, 0 evaluates the seperable kernel, 1 its 
          %  derivative and 2 its second derivative. Default is 0.
          %----------------------------------------------------------------
          % OUTPUT
          % vals  a obj.D by obj.D cell array containing an SepKernel object
          %       for each element of the Hessian matrix of the inputed
          %       SepKernel.
          %
          %----------------------------------------------------------------
          % EXAMPLES
          %% %% D = 2 
          % %% % General example
          % % generate a separable Gaussian kernel
          % gK = SepKernel( 2, [ 3, 6 ] );
          % % Define xvals vector
          % xvals = cell( [1 2] );
          % xvals{1} = -23:0.5:15;
          % xvals{2} = -5:0.25:3;
          % % evaluate the kernel and plot
          % evals = eval( gK, xvals, 0 );
          % figure, clf,
          % imagesc(evals)
          % %% % Default value for dx numeric
          % % evaluate the kernel and plot
          % evals = eval( gK, 1, 0 );
          % figure, clf,
          % imagesc( evals )
          % title( "dx=1" )
          % % evaluate the kernel and plot
          % figure, clf,
          % evals = eval( gK, 0.5, 0 );
          % imagesc( evals )
          % title( "dx=0.5" )
          % % evaluate the kernel and plot
          % figure, clf,
          % evals = eval( gK, 0.25, 0 );
          % imagesc( evals )
          % title( "dx=0.25" )
          %
          %----------------------------------------------------------------
          % Author: Fabian Telschow
          %----------------------------------------------------------------
          %%% Check mandatory input
          
          if iscell( xvals )
              if length( xvals ) == 1
                  tmp = xvals{1};
                  xvals = cell( [ 1 obj.D ] );
                  for d = 1:obj.D
                      xvals{d} = tmp;
                  end
              elseif length( xvals ) ~= obj.D
                  error( "Please, provide xvalues for each dimension." )
              end
          elseif numel( xvals ) == 1
              dx = xvals;
              xvals = cell( [ 1 obj.D ] );
              for d = 1:obj.D
                  if all( size( obj.truncation ) == [ 2 obj.D ] )
                      xvals{d} = -obj.truncation(1,d):dx:obj.truncation(2,d);
                  elseif all( size( obj.truncation ) == [ 1 obj.D ] )
                      xvals{d} = -obj.truncation(d):dx:obj.truncation(d);
                  else
                      error("Your truncation field in your SepKernel is incorrectly specified.")
                  end
              end
          else
              tmp = xvals;
              xvals = cell( [ 1 obj.D ] );
              for d = 1:obj.D
                  xvals{d} = tmp;
              end
          end
          
          %%% Main function
          switch obj.D
              case 1
                  locs1 = xvals{1};
                  locs{1} = locs1(:);
              case 2
                  % Get the meshgrid
                  [ locs1, locs2 ] = meshgrid( xvals{1}, xvals{2} );

                  % Make the locations into a cell structure
                  locs    = cell( [ 1 obj.D ] );
                  locs{1} = locs1(:);
                  locs{2} = locs2(:);
                  
                  % Get the size of the meshgrid
                  slocs = size( locs1 );
              case 3
                  % Get the meshgrid
                  [ locs1, locs2, locs3 ] = meshgrid( xvals{1},...
                                                      xvals{2}, xvals{3} );

                  % Make the locations into a cell structure
                  locs    = cell( [ 1 obj.D ] );
                  locs{1} = locs1(:);
                  locs{2} = locs2(:);
                  locs{3} = locs3(:);
                  
                  % Get the size of the meshgrid
                  slocs = size( locs1 );
          end
          
          % Evaluate the appropriate kernel on the grid
          switch derivtype
              case 0
                  % Evaluate the kernel on the meshgrid
                  vals = ones( [ prod( slocs ) 1 ] );
                  for d = 1:obj.D
                      vals = vals .* obj.kernel{d}( locs{d}(:) );
                  end

                  % Shape output back toslocs output size
                  if obj.D >1
                    vals = reshape( vals, slocs );
                  end
              case 1
              case 2
          end
      end
           
   end
end