classdef Field
   % Field is a class to make the use of multidimensional arrays, which are
   % considered to represent a field over a domain easy.
   %   The class is an object keeping track of the field and further
   %   important properties of the field.
   properties
      field double      % an array representing a field
      mask  { mustBeLogical( mask ) } = true % a logical array indicating which elements of the domain are inside the domain of the field
      xvals cell % a cell of size 1xlength(size(mask)) containing the coordinates of the voxels for each dimension
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
       %%% Validate compatibility of the properties
       % Note this solution is kind of a work around since validation
       % functions in matlab seem to work, badly. Hence we simply redefine
       % the set.property function.       
       % Change set.mask function
       function obj = set.mask( obj, val )
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
               sVal = size( val );
               l = length( sM );
               if ~all( sVal( 1:l ) == sM )
                  error( 'mask needs to have a compatible dimension with field.' )
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
                     elseif ~size( val{d}, 1 )
                         error( 'Entries of xvals needs to be a 1xL vector.' )
                     end
                end       
                % Assign the value
                obj.xvals = val;
             
           elseif mask_set % mask field provided
               sM = size( obj.mask );
               D  = length( val );
               l  = length( sM );
               if D == l
                   for d = 1:D
                         if ~isnumeric( val{d} )
                             error( 'xvals entries needs to be a numeric vector.' )
                         elseif ~size( val{d}, 1 )
                             error( 'Entries of xvals needs to be a 1xL vector.' )
                         elseif length( val{d} ) ~= sM(d)
                                 error('xvals needs to be compatible with mask.')
                         end
                   end
                   % Assign the value
                   obj.xvals = val;   
               end

           else
               D  = length( val );
               sF = size( obj.field );
               l  = length( sF );
               if D <= l
                   for d = 1:D
                       if ~isnumeric( val{d} )
                           error( 'xvals entries needs to be a numeric vector.' )
                       elseif ~size( val{d}, 1 )
                           error( 'Entries of xvals needs to be a 1xL vector.' )
                       elseif length( val{d} ) ~= sF(d)
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
        
       %%% Fill the dependent properties
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
       
       %%% Functions for class Field
       % Function for masking data
       obj = Mask( obj, mask )
       
   end
end

       function [] = mustBeLogical( F )
         if ~all( islogical( F(:) ) )
             error( [ 'Value assigned to mask property needs to be a ',...
                      'logical array'  ] )
         end
       end