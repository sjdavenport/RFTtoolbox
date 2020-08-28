classdef ConvFieldParams
   % ConvFieldParams is a class storing the parameters for generation of
   % a convolution field.
   properties
      kernel  SepKernel % an object of class SepKernel, which was used to generate the convolution field. 
      resadd(1,1)  { mustBeNonnegative, mustBeInteger } % an integer denoting how many voxels are introduced between two existing voxels.
      enlarge(1,1) { mustBeNonnegative, mustBeInteger } % an integer denoting by how many voxels in the high resolution the original boundary of the mask is dilated. 
      lat_masked(1,1) % a logical denoting whether the convolution field should first apply a mask.
   end
   methods
       %% Fill the dependent properties
       %-------------------------------------------------------------------
       
       %% Functions for class ConvField
       %-------------------------------------------------------------------
       % Change set.mask function
       function obj = set.lat_masked( obj, val )
           if ~islogical( val )
               error( "mask must be an logical array." )
           else
               obj.lat_masked = val;
           end
       end
       
       %% Basic constructor
       %-------------------------------------------------------------------
       function obj = ConvFieldParams( varargin )
          % ConvFieldParams( varargin ) is a basic constructor for a 
          % ConvFieldParams class object.
          % Default value for resadd is 1, for enlarge is ceil(resadd/2) and
          % for lat_masked is true.
          %----------------------------------------------------------------
          % ARGUMENTS
          % varargin
          %  If 1: Input must be a SepKernel object or a numeric vector
          %  If 2: 
          %     - first input is a SepKernel object or a numeric vector
          %     - second is an integer representing resadd 
          %  If 3:
          %     - first input is a SepKernel object or a numeric vector
          %     - second is an integer representing resadd 
          %     - third is an integer representing enlarge 
          %  If 4:
          %     - first input is a SepKernel object or a numeric vector
          %     - second is an integer representing resadd 
          %     - third is an integer representing enlarge 
          %     - forth is an logical representing lat_masked 
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
          if isnumeric( varargin{1} )
             obj.kernel = SepKernel( length( varargin{1} ), varargin{1} );
          else
             obj.kernel = varargin{1};
          end
          
          if nargin >= 2
              obj.resadd = varargin{2};
          else
              obj.resadd = 1;
          end
          
          if nargin >= 3
              obj.enlarge = varargin{3};
          else
              if mod( obj.resadd, 2 ) == 1
                obj.enlarge = ceil( obj.resadd / 2 );
              else
                obj.enlarge = 0;
              end
          end

          if nargin == 4
              obj.lat_masked = varargin{4};
          else
              obj.lat_masked = true;
          end
      end
   end
end