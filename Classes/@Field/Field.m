classdef Field
   % Field is a class to make the use of multidimensional arrays, which are
   % considered to represent a field over a domain easy.
   %   The class is an object keeping track of the field and further
   %   important properties of the field.
   properties
      field double % an array representing a field
      mask         % a logical array indicating which elements of the domain are inside the domain of the field
      xvals cell   % a cell of size 1xlength(size(mask)) containing the coordinates of the voxels for each dimension
   end
   properties ( Dependent ) % @Sam: Are these all properties we care about
                            % derived from the domain and the field?
                            % I think it would be useful to add the xvec
                            % here as well! What do you think?
      sizeField  % the size of the Field
      Dim        % the dimension of the domain 
      sizeDomain % the size of the domain of the field
      fiberDim   % the dimension of the fiber above a point
      sizeFiber  % the size of the fiber above a point
      masked     % 1 if the field is masked, zero else
   end
   properties ( Dependent, Access = private ) 
      dim_xvals  % dimension obtained from the provided xvals vector
   end
   methods
       %% Validate compatibility of the properties
       % Note this solution is kind of a work around since validation
       % functions in matlab seem to work, badly. Hence we simply redefine
       % the set.property function.
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
               if ~all( sM == obj.dim_xvals )
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
           
           if ~mask_set  && isempty( tmp.xvals ) % neither mask nor xvals specified
             % Assign the value
             obj.field = val;
           elseif mask_set % mask already specified
               sM = size( obj.mask );
               sF = size( val );
               l  = length( sM );
               ll = length( sF );
               if ll < l
                   error( "field must have at least as many dimensions as mask." )
               else
                   if ~all( sF( 1:l ) == sM )
                      error( 'field must have a compatible dimension with mask.' )
                   end
               end
               % Assign the value
               obj.field = val;
           else
               sF = size( val ); 
               if ~all( sF( 1:D ) == obj.dim_xvals )
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
                D  = length( val );
                for d = 1:D
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
               D  = length( val(:) );
               l  = length( sM );
               if D == l
                   for d = 1:D
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
               D  = length( val );
               sF = size( obj.field );
               l  = length( sF );
               if D <= l
                   for d = 1:D
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
       % Fill the sizeField field
       function value = get.sizeField( obj )
         value = size( obj.field );
       end
       
       % Fill the sizeDomain field
       function value = get.sizeDomain( obj )
         value = size( obj.mask );
       end
       
       % Fill the D field
       function D = get.Dim( obj )
         D = length( obj.sizeDomain );
         if any( obj.sizeDomain == 1) && D == 2
             D = 1;
         end
       end
       
       % Fill the fiberDim field
       function D = get.fiberDim( obj )
         D = length( obj.sizeField ) - length( obj.sizeDomain );
         if D == 0
             D = 1;
         end
       end
       
       % Fill the sizeFiber field
       function sF = get.sizeFiber( obj )
           if obj.fiberDim == 1
               sF = obj.sizeField( end );
               
               if length( obj.sizeField ) == obj.Dim
                   sF = 1;
               end
               
           else
               sF = obj.sizeField( (obj.Dim+1):end );
           end
       end
       
       % Fill the masked field
       function masked = get.masked( obj )
           tmp = obj.field( ~obj.mask(:) );
           
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

       %% Basic constructor
       %-------------------------------------------------------------------
       function obj = Field( varargin )
          % FIELD( varargin ) is a basic constructor for a Field class
          % object.
          %----------------------------------------------------------------
          % ARGUMENTS
          % varargin
          %  If 1:  
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
          if nargin == 0
             obj.field = [];
             obj.mask = true;
             obj.xvals = {};
          elseif nargin == 1
              if isnumeric( varargin{1} )
                % Check that varargin{1} is a vector
                Dim = varargin{1};
                sDim = size( Dim );
                if ~( length( sDim ) == 2 && sDim(1) == 1 )
                    error( "Input must be a 1xD vector." )
                end
                
                % Set mask to be all true
                obj.mask = true( Dim );
             
                % Set xvals to be equidistant in all dimensions with 
                % voxelsize 1
                xvals = cell( [ 1 obj.Dim ] );
                for d = 1:obj.Dim
                     xvals{d} = 1:obj.sizeDomain(d);
                end
                obj.xvals = xvals;
                
              elseif iscell( varargin{1} )
                  obj.mask  = true;
                  obj.xvals = varargin{1};
                  obj.mask  = true( obj.dim_xvals );
                  
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
       
       %%% Functions for class Field
       % Function for masking data
       obj = Mask( obj, mask )
       
   end
end