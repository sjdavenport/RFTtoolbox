classdef Field
   % Field is a class to make the use of multidimensional arrays, which are
   % considered to represent a field over a domain easy.
   %   The class is an object keeping track of the field and further
   %   important properties of the field like dimensions of the domain and 
   %   the fiber as well as their size. Moreover, it ensures that all these
   %   quantities are consistent.
   %   A complete field object consists of the following fields:
   %   - field   a T_1 x ... x T_D x F_1 x ... x F_K numeric array
   %             containing the values of the field. The T_1 x ... x T_D
   %             part is assumed to be the domain of the field and the
   %             part F_1 x ... x F_K is considered the fiber. A field
   %             of dimension T_1 x ... x T_D x {1} is a scalar field.
   %             Multiple observations of a scalar field are
   %             represented as a T_1 x ... x T_D x F_1 field, where
   %             F_1 denotes the sample size.
   %   - mask    a T_1 x ... x T_D logical array. Denoting a
   %             restriction of the field to a reduced domain given by
   %             the true values in mask. The dependend field 'masked'
   %             indicates whether the 'field' field has only all NaN,
   %             0 or -Inf for the values -Inf. Transforming a
   %             non-masked to a masked Field object can be achieved
   %             by the function Masked().
   %   - xvals   a 1 x D cell array which contains in the d-th entry
   %             the grid coordinates for the field. Note that a field
   %             object currently assumes that it is defined over a
   %             grid, i.e. it is defined on the cartesian product 
   %             xvals{1} x ... x xvals{D}.
   properties
      field  {mustBeNumeric} % an array representing a field
      mask         % a logical array indicating which elements of the domain are inside the domain of the field
      xvals cell   % a cell of size 1 x length(size(mask)) containing the coordinates of the voxels for each dimension
   end
   properties ( Dependent, Access = public ) 
      fieldsize  % the size of the Field
      D        % the dimension of the domain 
      masksize % the size of the domain of the field
      fiberD   % the dimension of the fiber above a point
      fibersize  % the size of the fiber above a point
      masked     % 1 if the field is masked, zero else
   end
   properties ( Dependent, Access = private ) 
      dim_xvals  % dimension obtained from the provided xvals vector
      complete   % indicates whether 'field', 'mask' and 'xvals' fields are
                 % all set to compatible values.
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
           if isempty( obj.field ) && isempty( obj.xvals )
             % Assign the value
             obj.mask = val;               
           elseif ~isempty( obj.field )
               sF = size( obj.field );
               sM = size( val );
               l = length( sM );
               if ~all( sM == sF( 1:l ) )
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
          % input the field 'field' is always empty and needs to be
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
          % % create an isotropic Gaussian kernel object
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
                  obj.mask  = true( obj.dim_xvals );
                  
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
       
       %% Functions for class Field
       %-------------------------------------------------------------------

       % Function for masking data
       obj = Mask( obj, mask )
       
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
            out.Kernel = Kernel;
            out.resadd = resadd;
            out.enlarge = enlarge;
            out.field = obj.field;
            out.mask = obj.mask;
            out.xvals = obj.xvals;            
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
       
       % Redefine plot
       out = plot( field, slice )
           
       % Redefine imagesc
       out = imagesc( field, slice, subj )
       
   end
end