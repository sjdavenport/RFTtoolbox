classdef VoxManifold
   % VoxManifold is a class, which summarizes geometrical information of a
   % manifold defined by voxels.
   %    It contains the mask, i.e., the points belonging to the manifold.
   %    Note that the code is ambigious in the sense that the boundary
   %    voxels of the mask are considered to be on the boundary.
   properties
      mask                 % a logical array indicating which elements define a voxel
      xvals        cell    % the size of the voxel in the different dimensions
      g            Field   % an object of class Field representing the Riemannian metric
      Gamma        Field   % an object containing the Christoffel symbols with respect to g
   end
   properties ( Dependent ) 
      size     % the size of the Field
      D        % the dimension of the voxel manifold
      masked   % 1 if the field objects are masked, zero else
   end
   properties ( Dependent, Access = private ) 
      dim_xvals  % dimension obtained from the provided xvals vector
   end
   methods
       %% Validate compatibility of the properties
       % Note this solution is kind of a work around since validation
       % functions in matlab seem to work, badly. Hence we simply redefine
       % the set.property functions to ensure compatibility between the
       % fields of a Field object.
       %-------------------------------------------------------------------
       % Change set.mask function
       function obj = set.mask( obj, val )
           if ~islogical( val )
               error( "mask must be an logical array." )
           end
           if isempty( obj.g.field ) && isempty( obj.xvals )
             % Assign the value
             obj.mask = val;               
           elseif ~isempty( obj.g.field )
               % Create a Field object from the input
               smask = size( val );
               if ~all( smask == obj.g.masksize ) && ~isempty( obj.g.field )
                  error('mask needs to have a compatible dimensions with g.')
               end
               % Assign the value
               obj.mask = val;
           else
               sM = size( val ); 
               if ~all( sM(1:obj.D) == obj.dim_xvals )
                  error( 'mask needs to have a compatible dimension with xvals.' )
               end
               
               % Assign the value
               obj.mask = val;               
           end
       end

       % Change set.g function
       function obj = set.g( obj, val )
           % 1 if mask is already filled
           mask_set = ~( length( obj.mask(:) ) == 1 ) && ~isempty( obj.mask );
           
           if ~mask_set  && isempty( obj.xvals ) % neither mask nor xvals specified
               if all( val.fibersize == val.D*ones( [ 1  2 ] ) )
                    % Assign the value
                    obj.g = val;
               else
                   error( "The fiber of the Riemannian metric needs to be D x D." )
               end
           elseif mask_set % mask already specified
               % Create a Field object from the input
               smask = size( obj.mask );
               if ~all( smask == val.masksize ) && ~isempty( obj.g.field )
                  error('mask needs to have a compatible dimensions with g.')
               end
               if all( val.fibersize == obj.D * [ 1 1 ] )
                    % Assign the value
                    obj.g = val;
               else
                   error( "The fiber of the Riemannian metric needs to be D x D." )
               end

           else
               sF = size( val ); 
               if ~all( sF( 1:obj.D ) == obj.dim_xvals )
                  error( 'fields needs to have a compatible dimension with xvals.' )
               end
               
               if all( val.fibersize == obj.D * [ 1 1 ] )
                    % Assign the value
                    obj.g = val;
               else
                   error( "The fiber of the Riemannian metric needs to be D x D." )
               end
           end
       end
       
       % Change set.xvals function
       function obj = set.xvals( obj, val )
           mask_set = ~( length( obj.mask(:) ) == 1 && obj.mask );
           if ~mask_set && isempty( obj.g.field ) % other fields not provided
                % Catch wrong xvals input  
                DD  = length( val );
                for d = 1:DD
                     if ~isnumeric( val{d} )
                         error( 'xvals entries needs to be a numeric vector.' )
                     end
                     if size( val{d}, 1 ) ~= 1
                         error( 'Entries of xvals needs to be a 1xL vector.' )
                     end
                end       
                % Assign the value
                obj.xvals = val;
             
           elseif mask_set % mask field provided
               sM = size( obj.mask );
               DD  = length( val(:) );

               if DD == obj.D
                   for d = 1:DD
                         if ~isnumeric( val{d} )
                             error( 'xvals entries needs to be a numeric vector.' )
                         end
                         if size( val{d}, 1 ) ~= 1
                             error( 'Entries of xvals needs to be a 1xL vector.' )
                         end
                         if length( val{d} ) ~= sM(d)
                                 error( 'xvals needs to be compatible with mask.' )
                         end
                   end
                   % Assign the value
                   obj.xvals = val;
               else
                   error( 'xvals needs to have a cell for each dimension in mask.' )
               end

           else
               DD  = length( val );
               g = obj.g;
               if DD == g.D
                   for d = 1:DD
                       if ~isnumeric( val{d} )
                           error( 'xvals entries needs to be a numeric vector.' )
                       end
                       if size( val{d}, 1 ) ~= 1
                           error( 'Entries of xvals needs to be a 1xL vector.' )
                       end
                       if length( val{d} ) ~= sF(d)
                             error('xvals needs to be compatible with mask.')
                       end
                   end
                   % Assign the value
                   obj.xvals = val;                    
               else
                    error('length(xvals) must be equal to the dimension of g.')
               end
           end

       end
        
       %% Fill the dependent properties
       %-------------------------------------------------------------------
       % Fill the D field
       function value = get.D( obj )
         % If the Riemannian metric property exist
         if ~isempty( obj.g.xvals )
            value = obj.g.D;

         % If the mask property exist
         elseif length( obj.mask(:) ) ~= 1
            sM     = size( obj.mask );
            value  = length( sM );
            if value == 2 && sM(2) == 1
                value = 1;
            end

         % If the xvals property exist
         elseif ~isempty( obj.xvals )
            value = length( obj.xvals );
         end
       end
       
       % Fill the size field
       function value = get.size( obj )
            value = size( obj.mask );
       end
       
       % Fill the masked field
       function masked = get.masked( obj )
           % Check whether the Field object obj.g is masked using the mask
           % property
           s = prod( obj.g.fibersize );
           tmp = obj.g.field( repmat( ~obj.mask(:), [ s 1 ] ) );
           
           masked = 0;
           if all( tmp == 0 ) || all( isnan( tmp ) ) || all( tmp == -Inf )
               masked = 1;
           end
       end

       % Fill the dim_xvals field
       function dim_xvals = get.dim_xvals( obj )
           DD  = length( obj.xvals );
           dim_xvals = zeros( [ 1 DD ] );
           for d = 1:DD
               dim_xvals(d) = length( obj.xvals{d} );
           end
       end

       %% Basic constructor
       %-------------------------------------------------------------------
       function obj = VoxManifold( varargin )
          % VoxManifold( varargin ) is a basic constructor for an object of
          % class VoxManifold.
          % Default values for mask and xvals are given by g.
          %----------------------------------------------------------------
          % ARGUMENTS
          % varargin
          %  If 0: Generates empty field object.
          %  If 1: possible inputs are
          %     - an complete object of class Field defining the Riemannian
          %       metric. Mask, xvals and masked property are inherited
          %       from g.
          %     - 1 x D vector defining the size of the mask. 'mask' and 
          %       'xvals' are set to the default value.
          %     - 1 x D cell array containing the xvals. 'mask' is default.
          %     - T_1 x ... x T_D logical array containing the mask.
          %       'xvals' is set to default.
          %  If 2: possible inputs are
          %     - T_1 x ... x T_D x F_1 x ... x F_K numerical array, which
          %         defines the values of the field.
          %       T_1 x ... x T_D logical array containing the mask.
          %       'xvals' is set to default.
          %     - T_1 x ... x T_D x F_1 x ... x F_K numerical array, which
          %         defines the values of the field.
          %       1 x D cell array containing the xvals.
          %       'mask' is default.
          %     - T_1 x ... x T_D x F_1 x ... x F_K numerical array, which
          %         defines the values of the field.
          %       an integer D, which specifies 'D'
          %       'mask' and 'xvals' are default.
          %  If 3: possible inputs are
          %     - T_1 x ... x T_D x F_1 x ... x F_K numerical array, which
          %         defines the values of the field.
          %       T_1 x ... x T_D logical array containing the mask.
          %       1 x D cell array containing the xvals.
          %----------------------------------------------------------------
          % OUTPUT
          % obj  an object of class Field where the fields are
          %      initialised by the input or defaults.
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
              if isa( varargin{1}, 'Field' ) % Input is the Riemannian metric
                if iscompletefield( varargin{1} )
                    obj.g = varargin{1};
                    obj.mask = obj.g.mask;
                    obj.xvals = obj.g.xvals;
                else
                    error( "The input Field object needs to be complete, i.e. all entries must be filled." )
                end

              elseif isnumeric( varargin{1} ) % Input is vector
                  obj.g     = EuclideanMetric( Field( true( varargin{1} ) ) );
                  obj.mask  = true( varargin{1} );
                  % Default the Riemannian metric to be Euclidean
                  obj.g = EuclideanMetric( obj.g, obj.mask );
                  % Set xvals to be equidistant in all dimensions with 
                  % voxelsize 1
                  xvals = cell( [ 1 obj.D ] );
                  for d = 1:obj.D
                       xvals{d} = 1:varargin{1}(d);
                  end
                  obj.xvals = xvals;
                  
                  % Default Christoffelsymbols are zero
                  obj.Gamma = constfield( 0, [ obj.masksize obj.D * [ 1 1 1] ],...
                                                                obj.xvals );
                  obj.Gamma.mask = obj.mask;
                
              elseif iscell( varargin{1} ) % Input is xvals
                  obj.g       = EuclideanMetric( Field( varargin{1} ) );
                  obj.g.xvals = varargin{1};
                  obj.mask  = obj.g.mask;
                  obj.xvals = varargin{1};
                  obj.mask  = true( obj.dim_xvals );
                  obj.g     = EuclideanMetric( obj.g, obj.mask );
                  % Default Christoffelsymbols are zero
                  obj.Gamma = constfield( zeros( obj.D * [ 1 1 1 ] ), obj.size ,...
                                                                obj.xvals );
                  obj.Gamma.mask = obj.mask;
                  
              elseif islogical( varargin{1} ) % Input is mask
                  obj.g       = EuclideanMetric( Field( varargin{1} ) );
                  obj.mask  = varargin{1};
                  % Set xvals to be equidistant in all dimensions with 
                  % voxelsize 1
                  xvals = cell( [ 1 obj.D ] );
                  for d = 1:obj.D
                       xvals{d} = 1:obj.size(d);
                  end
                  obj.xvals = xvals;
                  % Default Christoffelsymbols are zero
                  obj.Gamma = constfield( zeros( obj.D * [ 1 1 1 ] ), obj.size ,...
                                                                obj.xvals );
                  obj.Gamma.mask = obj.mask;
                  
              else
                  error( strcat( "Input must be either a cell",...
                                 "containing the locations or a numerical vector" ) )
              end
              
          elseif nargin == 2
              if islogical( varargin{1} ) && iscell( varargin{2} ) % Input is mask and xvals
                  obj.g     = EuclideanMetric( Field( varargin{1} ) );
                  obj.mask  = varargin{1};
                  obj.xvals = varargin{2};
                  obj.g     = EuclideanMetric( obj.g, obj.mask );
                  % Default Christoffelsymbols are zero
                  obj.Gamma = constfield( zeros( obj.D * [ 1 1 1 ] ), obj.size ,...
                                                                obj.xvals );
                  obj.Gamma.mask = obj.mask;
                  
              elseif isa( varargin{1}, 'Field' ) && isnumeric( varargin{2} ) % Input is g and resadd
                  obj         = VoxManifold( varargin{1} );
                  % Default Christoffelsymbols are zero
                  obj.Gamma = constfield( zeros( obj.D * [ 1 1 1 ] ), obj.size ,...
                                                                obj.xvals );
                  obj.Gamma.mask = obj.mask;
              else
                  error( strcat( "VoxManifold constructor not available",...
                                 " for the chosen input. Please refer to",...
                                 " the description of the class." ) );
              end
          elseif nargin == 3
              if islogical( varargin{1} ) && iscell( varargin{2} )...
                                            && isnumeric( varargin{3} )
                  % Input is mask, xvals and resadd
                  obj.g     = EuclideanMetric( Field( varargin{1} ) );
                  obj.mask  = varargin{1};
                  obj.xvals = varargin{2};
                  % Default Christoffelsymbols are zero
                  obj.Gamma = constfield( zeros( obj.D * [ 1 1 1 ] ), obj.size ,...
                                                                obj.xvals );
                  obj.Gamma.mask = obj.mask;
                  
              elseif isa( varargin{1}, 'Field' ) && isnumeric( varargin{2} )...
                                                    && isnumeric( varargin{3} )
                  % Input is g, resadd and enlarge
                  obj         = VoxManifold( varargin{1} );
                  % Default Christoffelsymbols are zero
                  obj.Gamma = constfield( zeros( obj.D * [ 1 1 1 ] ), obj.size ,...
                                                                obj.xvals );
                  obj.Gamma.mask = obj.mask;
                  
              else
                  error( strcat( "VoxManifold constructor not available",...
                                 " for the chosen input. Please refer to",...
                                 " the description of the class." ) );
              end
          else
                error( "Too many input variables." )
          end
      end
       
       %% Functions for class VoxManifold
       %-------------------------------------------------------------------
       % Inner integral in the \partial_1 M part of LKC1
       angles = IntegralAlpha( voxmfd, mask, eucangle )
       
       % Estimation of Lipschitz-Killing curvatures
       [ L, L0 ] = LKC_est( voxmfd, version )
       
       % Inner product between vector fields
       F = InnerProd( voxmfd, v, w )
       
   end
end