classdef Field
   % Field is a class making the use of multidimensional arrays, which are
   % considered to represent a field (or for mathematicians a fiber bundle)
   % over a domain.
   %   The class is an object keeping track of the field and important
   %   properties thereof, for example, the dimensions of the domain
   %   and the fiber as well as their sizes. Moreover, it ensures that all
   %   these quantities are consistent.
   %   A complete Field object consists of the following properties:
   %   - field   a T_1 x ... x T_D x F_1 x ... x F_K numeric array
   %             containing the values of the field. The T_1 x ... x T_D
   %             part is assumed to be the domain of the field and the
   %             part F_1 x ... x F_K is considered the fiber. A field
   %             of dimension T_1 x ... x T_D x {1} is a scalar field.
   %             Multiple observations of a scalar field are
   %             represented as a T_1 x ... x T_D x F_1 field, where
   %             F_1 > 1 denotes the number of observations.
   %   - mask    a T_1 x ... x T_D logical array. Denoting a
   %             restriction of the field to a reduced domain given by
   %             the true values in mask. The dependend field 'masked'
   %             indicates whether the field property has only all NaN,
   %             0 or -Inf for the values outside the mask. Transforming a
   %             non-masked to a masked Field object can be achieved
   %             by the function Masked().
   %   - xvals   a 1 x D cell array which contains in the d-th entry
   %             the grid coordinates for the field. Note that a field
   %             object currently assumes that it is defined over a
   %             grid, i.e. it is defined on the cartesian product 
   %             xvals{1} x ... x xvals{D}.
   properties
      field  {mustBeNumeric} % an T_1 x ... x T_D x F_1 x ... x F_K numeric array representing the values of the field.
      mask         % a T_1 x ... x T_D logical array indicating which voxels belong to the domain of the field.
      xvals cell   % a 1 x D cell containing the d-th coordinate of each the voxels.
   end
   properties ( Dependent, Access = public ) 
      fieldsize  % the size of the Field, i.e., [ T_1, ..., T_D, F_1, ..., F_K ].
      D          % the dimension of the domain, i.e., D.
      masksize   % the size of the domain of the field, i.e., [ T_1, ..., T_D ].
      fiberD     % the dimension of the fiber above a point, i.e., K.
      fibersize  % the size of the fiber above a point, i.e., [ F_1, ..., F_K ].
      masked     % 1, if the field is masked, else 0.
   end
   properties ( Dependent, Access = private ) 
      dim_xvals  % dimension obtained from the provided xvals vector
      complete   % indicates whether 'field', 'mask' and 'xvals' properties
                 % are all set to compatible non-default values.
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
           % Logical input
           if ~islogical( val )
               error( "mask must be an logical array." )
           end
           
           % Catch wrong 1D mask input
           sM = size(val);
           if sM(1) == 1 && sM(2) ~= 1
               error( "a 1D mask must be a logical column vector." )
           end
           
           % Check compatibilitiy with field
           if isempty( obj.field ) && isempty( obj.xvals )
             % Assign the value
             obj.mask = val;               
           elseif ~isempty( obj.field )
               sF = size( obj.field );
               sM = size( val );
               sM = sM( 1:obj.D );
               
               if ~all( sM == sF( 1:obj.D ) )
                  error('mask needs to have a compatible dimension with field.')
               end
               % Assign the value
               obj.mask = val;
           else
               sM = size( val ); 
               if ~all( sM(1:obj.D) == obj.dim_xvals )
                  error( 'fields needs to have a compatible dimension with xvals.' )
               end
               
               % Assign the value
               obj.mask = val;               
           end
       end

       % Change set.field function
       function obj = set.field( obj, val )
           % 1 if mask is already filled
           mask_set = ~( length( obj.mask(:) ) == 1 );
           
           if ~mask_set  && isempty( obj.xvals ) % neither mask nor xvals specified
             % Assign the value
             obj.field = val;
           elseif mask_set % mask already specified
               sM = size( obj.mask );
               sF = size( val );
               l  = length( sM );
               ll = length( sF );
               if l > ll
                   error( "field must have at least as many dimensions as mask." )
               else
                   if ~all( sF( 1:obj.D ) == sM( 1:obj.D ) )
                      error( 'field must have a compatible dimension with mask.' )
                   end
               end
               % Assign the value
               obj.field = val;
           else
               sF = size( val ); 
               if ~all( sF( 1:obj.D ) == obj.dim_xvals )
                  error( 'fields needs to have a compatible dimension with xvals.' )
               end
               
               % Assign the value
               obj.field = val;
           end
       end
       
       % Change set.xvals function
       function obj = set.xvals( obj, val )
           mask_set = ~( length( obj.mask(:) ) == 1 && obj.mask );
           if ~mask_set && isempty( obj.field ) % other fields not provided
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
               sF = size( obj.field );
               l  = length( sF );
               if DD <= l
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
                    error('length(xvals) must be smaller or equal to the dimensions in field.')
               end
           end

       end
        
       %% Fill the dependent properties
       %-------------------------------------------------------------------
       % Fill the fieldsize field
       function value = get.fieldsize( obj )
         value = size( obj.field );
       end
       
       % Fill the masksize field
       function value = get.masksize( obj )
         value = size( obj.mask );
       end
       
       % Fill the D field
       function D = get.D( obj )
          if ~isempty( obj.xvals )
              D = length( obj.xvals );
          else 
              D = length( obj.masksize );
              if any( obj.masksize == 1) && D == 2
                  D = 1;
              end
          end
       end
       
       % Fill the fiberD field
       function D = get.fiberD( obj )
         l = length( obj.fieldsize );
         if l == 2 && obj.fieldsize( 2 ) == 1
             l = 1;
         end
         D = l - obj.D;
         if D == 0
             D = 1;
         end
       end
       
       % Fill the fibersize field
       function sF = get.fibersize( obj )
           if obj.fiberD == 1
               sF = obj.fieldsize( end );
               
               if length( obj.fieldsize ) == obj.D
                   sF = 1;
               end
               
           else
               sF = obj.fieldsize( (obj.D+1):end );
           end
       end
       
       % Fill the masked field
       function masked = get.masked( obj )
           s = prod( obj.fibersize );
           tmp = obj.field( repmat( ~obj.mask(:), [ s 1 ] ) );
           
           masked = 0;
           if all( tmp == 0 ) || all( isnan( tmp ) ) || all( tmp == -Inf )
               masked = 1;
           end
       end

       % Fill the dim_xvals field
       function dim_xvals = get.dim_xvals( obj )
           D  = length( obj.xvals );
           dim_xvals = zeros( [ 1 D ] );
           for d = 1:D
               dim_xvals(d) = length( obj.xvals{d} );
           end
       end
       
       % Fill the complete field
       function value = get.complete( obj )
           if ~isempty( obj.field ) && ~isempty( obj.xvals )
                value = true;
           else
               value = false;
           end
       end

       %% Basic constructor
       %-------------------------------------------------------------------
       function obj = Field( varargin )
          % FIELD( varargin ) is a basic constructor for a Field class
          % object.
          % Default values for mask is always true( masksize ) and default
          % for xvals is { 1:masksize(1), ..., 1:masksize(D) }. If not an
          % input the property 'field' is always empty and needs to be
          % specified by the user.
          %----------------------------------------------------------------
          % ARGUMENTS
          % varargin
          %  If 0: Generates empty field object.
          %  If 1: possible inputs are
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
          % % create an isotropic Gaussianfsu kernel object
          % D = 3;
          % gaussK = SepKernel( D, 6 )
          %
          %----------------------------------------------------------------
          % Author: Fabian Telschow
          %----------------------------------------------------------------
          if nargin == 0
             obj.field = [];
             obj.mask = true;
             obj.xvals = {};
          elseif nargin == 1
              if isnumeric( varargin{1} ) % Input is masksize
                % Check that varargin{1} is a vector
                Dim = varargin{1};
                sDim = size( Dim );
                if ~( length( sDim ) == 2 && sDim(1) == 1 )
                    error( "Input must be a 1xD vector." )
                end
                
                % Set mask to be all true
                obj.mask = true( [ Dim 1 ] );
             
                % Set xvals to be equidistant in all dimensions with 
                % voxelsize 1
                xvals = cell( [ 1 obj.D ] );
                for d = 1:obj.D
                     xvals{d} = 1:obj.masksize(d);
                end
                obj.xvals = xvals;
                
              elseif iscell( varargin{1} ) % Input is xvals
                  obj.mask  = true;
                  obj.xvals = varargin{1};
                  obj.mask  = true( [ obj.dim_xvals, 1] );
                  
              elseif islogical( varargin{1} ) % Input is mask
                  obj.mask  = varargin{1};
                  % Set xvals to be equidistant in all dimensions with 
                  % voxelsize 1
                  xvals = cell( [ 1 obj.D ] );
                  for d = 1:obj.D
                       xvals{d} = 1:obj.masksize(d);
                  end
                  obj.xvals = xvals;
                  
              else
                  error( strcat( "Input must be either a cell",...
                                 "containing the locations or a numerical vector" ) )
              end
          elseif nargin == 2
              if isnumeric( varargin{1} ) % First input is a field
                  if isnumeric( varargin{2} ) % Second input is the dimension of domain
                      D = varargin{2};
                      if length( D(:) ) == 1
                          % Find the size of the domain
                          sF = size( varargin{1} );
                          % Initialize standard field object of size sF
                          obj = Field( sF(1:D) );
                          % Fill the field field
                          obj.field = varargin{1};
                      else
                          error( "Second input must be an integer, i.e. the dimension of the domain." )
                      end
                      
                  elseif islogical( varargin{2} ) % Second input is the mask
                      % Get field of mask dimension
                      obj = Field( size( varargin{2} ) );
                      obj.mask  = varargin{2};
                      obj.field = varargin{1};
                  elseif iscell( varargin{2} ) % Second input is xvals
                      % Get field of mask dimension
                      obj = Field( varargin{2} );
                      obj.field = varargin{1};
                  else
                      error( "Input type for second input not supported." )
                  end
              else
                  error( "Input type for first input not supported." )
              end
          elseif nargin == 3
              obj = Field( varargin{1:2} );
              obj.xvals = varargin{3};
          else
                error( "Too many input variables." )
          end
      end
       
       %% General functions for class Field
       %-------------------------------------------------------------------
       function varargout = subsref(obj, s)
            if strcmp(s(1).type, '()')
                 % Catch insufficient bracket information
                 if length( s.subs ) ~= obj.D && ...
                         length( s.subs ) ~= obj.D + obj.fiberD
                     error( "You need to index either only the mask or the whole field." )
                 end

                 % Fill the fiber subs with ':', if only the mask is
                 % subsetted
                 if length( s.subs ) == obj.D
                     for k = 1:obj.fiberD
                        s.subs{ obj.D + k } = ':';
                     end
                 end
                     
                 ss = s;
                 ss.subs = s.subs( 1:obj.D );
                 newf = Field();
                 
                 % Make sure that the new mask is squeezed and a column in
                 % case of 1D slice
                 newmask = squeeze( builtin( 'subsref', obj.mask, ss ) );
                 if length( size( newmask ) ) == 2
                     if size( newmask, 1 ) == 1
                         newmask = newmask';
                     end
                 end
                 newf.mask = newmask;
                 
                 % Get new xvals structure
                 newxvals = obj.xvals;
                 l = NaN * ones( [ 1 obj.D ] );
                 for d = 1:obj.D
                     newxvals{d} = newxvals{d}( ss.subs{d} );
                     l(d) = length( newxvals{d} );
                     if strcmp( ss.subs{d}, ':' )
                         newxvals{d} = newxvals{d}';
                         l(d) = 666;
                     end
                 end
                 
                 % Remove dimensions with a single point
                 newxvals = newxvals( l ~= 1 );
                 newf.xvals   = newxvals;
                 
                 ff = squeeze( builtin( 'subsref', obj.field, s ) );
                 if size( ff, 1 ) == 1 && length( size( ff ) ) == 2
                     ff = ff';
                 end
                 newf.field   = ff;

                 varargout{1} = newf;
            else
                 [varargout{1:nargout}] = builtin('subsref', obj, s);
            end
       end
       
       % Get the stepsizes in each direction 
       dx = get_dx( obj )
        
       % Function for masking data
       obj = Mask( obj, val, mask )
       
       % Function for cutting a Field class to a mask
       cobj = cut2mask( obj )
       
       % Function for obtaining the private complete field
       function val = iscomplete( obj )
           val = obj.complete;
       end
       
       % Function for creating an Euclidean metric object
       obj = EuclideanMetric( obj, mask )
       
       % Transform a Field class object into a ConvField object
       function out = Field2ConvField( obj, Kernel, resadd, enlarge )
            if ~exist( 'Kernel', 'var' ) || isa( Kernel, 'SepKernel' )
                error( "An object of class SepKernel need to be provided." );
            end
            if ~exist( 'resadd', 'var' )
                resadd = 0;
            end
            if ~exist( 'enlarge', 'var' )
                enlarge = 0;
            end
            
            % Create an empty ConvField object
            out = ConvField();
            out.Kernel  = Kernel;
            out.resadd  = resadd;
            out.enlarge = enlarge;
            out.field   = obj.field;
            out.mask    = obj.mask;
            out.xvals   = obj.xvals;            
       end

       % Function for checking whether two Field objects are compatible
       function val = iscompatible( obj1, obj2 )
            val = false;
            if obj1.D == obj2.D
                if all( obj1.masksize == obj2.masksize )
                    val = true;
                    for d = 1:obj1.D
                        if ~all( obj1.xvals{ d } == obj2.xvals{ d } )
                            val = false;
                        end
                    end
                end
            end
       end
       
       %% Functions for Geometry
       %-------------------------------------------------------------------
       
       % Function converting two fields representing the field and its
       % derivative into a VoxelManifold
       voxmfd = Field2VoxManifold( field, dfield, masked )
       
       %% Plot functions for class Field
       %-------------------------------------------------------------------
       % Redefine plot
       out = plot( varargin )
           
       % Redefine imagesc
       out = imagesc( field, ntics )

       %% Stats Functions for class Field
       %-------------------------------------------------------------------            
       % Redefine sum()
       out = sum( varargin )
       
       % Redefine mean()
       out = mean( varargin )
       
       % Redefine std()
       out = std( varargin )
       
       % Redefine var()
       out = var( varargin )
       
       %% Algebra Functions for class Field
       %-------------------------------------------------------------------     
       % Redefine matrix multiplication
       out = mtimes( obj, v )
       
       % Redefine .* operator
       function out = times( obj1, obj2 )
           if isa( obj1, 'Field' )
               A = obj1.field;
               out = obj1;
           else
               A = obj1;
           end
           
           if isa( obj2, 'Field' )
               B = obj2.field;
               out = obj2;
           else
               B = obj2;
           end
           
           out.field = A .* B;
       end
       
       % Redefine ./ operator
       function out = rdivide( obj1, obj2 )
           if isa( obj1, 'Field' )
               A = obj1.field;
               out = obj1;
           else
               A = obj1;
           end
           
           if isa( obj2, 'Field' )
               B = obj2.field;
               out = obj2;
           else
               B = obj2;
           end
           
           out.field = A ./ B;
       end
       
       % Redefine .\ operator
       function out = ldivide( obj1, obj2 )
           if isa( obj1, 'Field' )
               A = obj1.field;
               out = obj1;
           else
               A = obj1;
           end
           
           if isa( obj2, 'Field' )
               B = obj2.field;
               out = obj2;
           else
               B = obj2;
           end
           
           out.field = A .\ B;
       end
       
       % Redefine .^ operator
       function out = power( obj1, obj2 )
           out = obj1;
           
           if isa( obj2, 'Field' )
               B = obj2.field;
               out = obj2;
           else
               B = obj2;
           end
                      
           out.field = obj1.field.^B;
       end
       
       % Redefine addition
       function out = plus( obj1, obj2 )
           if isa( obj1, 'Field' )
               A = obj1.field;
               out = obj1;
           else
               A = obj1;
           end
           
           if isa( obj2, 'Field' )
               B = obj2.field;
               out = obj2;
           else
               B = obj2;
           end
           
           out.field = A + B;
       end
       
       % Redefine subtraction
       function out = minus( obj1, obj2 )
           if isa( obj1, 'Field' )
               A = obj1.field;
               out = obj1;
           else
               A = obj1;
           end
           
           if isa( obj2, 'Field' )
               B = obj2.field;
               out = obj2;
           else
               B = obj2;
           end
           
           out.field = A - B;
       end

       % Redefine cross
       function obj = cross( obj1, obj2 )
           obj = obj1;
           obj.field = cross( obj1.field, obj2.field );
       end
       
       % Redefine sqrt
       function obj = sqrt( obj )
           obj.field = sqrt( obj.field );
       end
       
       % Redefine transpose
       function obj = transpose( obj )
           fD = obj.fiberD;
           if fD > 2
               error( "Transpose is only defined for a matrix or vector" )
           end
           
           if fD == 1
               obj.field = reshape( obj.field,...
                                    [ obj.masksize( 1:obj.D ), 1, obj.fibersize ] );
           elseif fD == 2 && obj.fibersize(1) == 1
               obj = squeeze( obj );
           else
           end
               
       end
       
       % Define inverse of symmetric 2x2 or 3x3 matrix
       function obj = invsym( obj )
           obj.field = invsym( obj.field );
       end
       
       % Define acos
       function obj = acos( obj )
           obj.field = acos( obj.field );
       end
       
       % Redefine equal
       function out = eq( obj1, obj2 )
           out = [ false, false, false ];
           
           if obj1.D == obj2.D && obj1.fiberD == obj2.fiberD
               if all( obj1.masksize == obj2.masksize ) ...
                        &&  all( obj1.fibersize == obj2.fibersize )
                   if all( obj1.field == obj2.field )
                       out(1) = true;
                   end

                   if all( obj1.mask == obj2.mask )
                       out(2) = true;
                   end

                   for d = 1:obj1.D
                        if all( obj1.xvals{d} == obj2.xvals{d} )
                            out(3) = true;
                        end 
                   end
               end
           end
       end
       
       % Define a subfield function
       function out = Subfield( obj, S )
            S1  = S(1:obj.D);
            out = Field( obj.mask( S1{:} ) );
            out.field = obj.field( S{:} );
            for d = 1:out.D
                val = obj.xvals{d}(S{d});
                if ischar( S{d} )
                    val = val';
                end
                out.xvals{d} = val;
            end
       end
       
       % Redefine a squeeze function
       function out = squeeze( obj )
            out       = Field( squeeze( obj.mask ) );
            out.field = squeeze( obj.field );
            ind       = false( [ 1 obj.D ] );
            for d = 1:obj.D
                if length( obj.xvals{d} ) > 1
                    ind( d ) = true;
                end
            end
            out.xvals = obj.xvals( ind );
       end
       
       % Collapse: a function collapsing either the fiber or the domain to
       % one dimension
       function out = collapse( obj, part )
            if ~exist( 'part', 'var' )
                part = "domain";
            end
            
            if strcmp( part, 'domain' )
                out = Field( obj.mask(:) );
                out.field = reshape( obj.field,...
                                     [ prod( obj.masksize ), obj.fibersize ] );
            elseif strcmp( part, 'fiber' )
                out = obj;
                out.field = reshape( obj.field,...
                                     [ obj.masksize, prod( obj.fibersize ) ] );
            else
                error( "You need to choose a valid part to collapse. Either domain or fiber." )
            end
       end
       
   end
end