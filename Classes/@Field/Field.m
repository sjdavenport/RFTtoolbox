classdef Field
   % Field is a class to make the use of multidimensional arrays, which are
   % considered to represent a field over a domain easy.
   %   The class is an object keeping track of the field and further
   %   properties important for the field.
   properties
      field       % a 1 by D cell array containing the function handles for the kernel in each direction
      mask        % a 1 by D array containing the value for truncation for the kernels
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
   methods
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
       
       % Function for masking data
       obj = Mask( obj, mask )
       
   end
end