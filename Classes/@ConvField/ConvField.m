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
   %                We recommend to only consider enlarge = 0, i.e.
   %                interpreting the boundary voxels of the original mask as
   %                being on the boundary of the domain of the convolution
   %                field. Another canonical choice is enlarge = ceil( resadd
   %                / 2 ), which means that each voxel of the original mask is
   %                considered to be in the center of a cube, which defines
   %                the manifold. In this case resadd needs to be odd in order
   %                to ensure that the boundary voxels of the high resolution
   %                mask are samples from the boundary.
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
       
       %% Functions for class Field
       %-------------------------------------------------------------------

       % Function converting a ConvField of derivetype 0 and a ConvField
       % of derivetype 1 into a VoxelManifold
       voxmfd = ConvField2VoxManifold( cfield, dcfield, masked )
       
       % Function to compute the Riemannian metric induced by a convolution
       % field
       g = Riemmetric_est( cfield, dcfield )
       
       % Function to compute the LKC from a voxel manifold obtained from
       % a convolution field
       [ L, L0 ] = LKC_voxmfd_est( cfield, dcfield )
   end
end