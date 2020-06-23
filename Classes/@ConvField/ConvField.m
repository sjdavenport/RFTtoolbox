classdef ConvField < Field
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
      lat_data
      resadd(1,1)  { mustBeNonnegative, mustBeInteger }
      enlarge(1,1) { mustBeNonnegative, mustBeInteger }
      derivtype(1,1) { mustBeNonnegative, mustBeInteger }
      kernel  SepKernel % a cell of size 1 x length(size(mask)) containing the coordinates of the voxels for each dimension
   end
   methods
       %% Functions for class Field
       %-------------------------------------------------------------------

       % Function converting a ConvField of derivetype 0 and a ConvField
       % of derivetype 1 into a VoxelManifold
       obj = ConvFields2VoxManifold( obj1, obj2 )
       
       % Function to compute the Riemannian metric induced by a convolution
       % field
       Lambda = Lambda_est( cfield, dcfield )
       
       % Function to compute the LKC from a voxel manifold obtained from
       % a convolution field
       [ L, L0 ] = LKC_voxmfd_est( cfield, dcfield )
   end
end