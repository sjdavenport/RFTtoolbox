classdef ConvField < Field
   % ConvField is a class storing fields derived from convolving lattice
   % data with a smoothing kernel.
   %   It is a subclass of the Field class and should be generated from a
   %   Field class object using the function convfield.m.
   %   Additional to Field it tracks information on its generation, i.e.:  
   %   - resadd     an integer denoting the amount of voxels added inbetween 
   %                each original voxel. It increases the resolution of the field
   %                by evaluation of the convolution field on these voxels.     
   %   - enlarge    an integer denoting by how many voxels in high resolution
   %                the mask is dilated. Note that this allows for different
   %                interpretations of the domains of the convolution field.
   %                We recommend to use enlarge = ceil( resadd
   %                / 2 ), which means that each voxel of the original mask
   %                is considered to be in the center of a cube, which defines
   %                the underlying manifold. In this case resadd needs to be
   %                odd in order to ensure that the boundary voxels of the 
   %                high resolution mask are samples from the boundary.
   %                Another canonical choice is enlarge = 0, i.e.
   %                interpreting the boundary voxels of the original mask as
   %                being on the boundary of the domain of the convolution
   %                field. This has the advantage that resadd can be
   %                arbitrary.
   %   - derivtype  an integer denoting whether the original field is
   %                stored or derivatives thereof.
   %   - kernel     an object of class SepKernel, which was used to generate
   %                the convolution field.
   properties
      resadd(1,1)  { mustBeNonnegative, mustBeInteger } % an integer denoting how many voxels are introduced between two existing voxels.
      enlarge(1,1) { mustBeNonnegative, mustBeInteger } % an integer denoting by how many voxels in the high resolution the original boundary of the mask is dilated. 
      derivtype(1,1) { mustBeNonnegative, mustBeInteger } % an integer denoting whether the smoothed field (0), its first derivatives (1) or its second derivatives (2) are stored in this field property.
      kernel  SepKernel % an object of class SepKernel, which was used to generate the convolution field. 
   end
   properties ( Dependent, Access = public ) 
      origsize % a 1 x D vector containing the size of the image before the resolution has been increased.
   end
   methods
       %% Fill the dependent properties
       %-------------------------------------------------------------------
       % Fill the fieldsize field
       function origsize = get.origsize( obj )
           enl = obj.enlarge;
           res = obj.resadd;
           sm  = obj.masksize;
           % Compute low resolution size
           origsize = sm - ( sm - 1 ) * res;
           % Correct for the enlargement
           if obj.D ==1
               origsize = origsize - [ 2 * enl, 0 ];
           else
               origsize = origsize - 2 * enl;
           end
       end
       
       %% Functions for class ConvField
       %-------------------------------------------------------------------
       function varargout = subsref(obj, s)
            if strcmp(s(1).type, '()')
                 ss = s;
                 ss.subs = s.subs(1:obj.D);
                 newf = ConvField();
                 newf.mask = builtin( 'subsref', obj.mask, ss );
                 newf.resadd = obj.resadd;
                 newf.xvals = obj.xvals;
                 newf.enlarge = obj.enlarge;
                 newf.field = builtin( 'subsref', obj.field, s );
                 varargout{1} = newf;
            else
                 [varargout{1:nargout}] = builtin('subsref', obj, s);
            end
       end
        
       % Function for checking whether two ConvField objects are compatible
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
                    if obj1.resadd ~= obj2.resadd || obj1.enlarge ~= obj2.enlarge
                        val = false;
                    end
                end
            end
       end       

   end
end