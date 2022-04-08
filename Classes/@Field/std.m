function out = std( varargin )
% STD( varargin ) does take the mean of specified fiber dimensions.
%--------------------------------------------------------------------------
% ARGUMENTS
% varargin
%  If 1: inputs are
%     - an object of class Field.
%  If 2: inputs are
%     - object of Class Field.
%     - positive integer denoting at which dimension in the fiber, the mean 
%       is taken.
%--------------------------------------------------------------------------
% OUTPUT
% obj  an object of class Field where the standard deviation of the
%      specified fiber is computed. If nargin==1, then the last
%      fiberdimension is taken the standard deviation of.
%
%--------------------------------------------------------------------------
% EXAMPLES
% % create a standard kernel object
% field = wfield( [3 3], [12 12 ] )
% % sum of last fiberdimension
% f1 = std( field )
% % sum over last  last fiberdimension different way
% f2 = std( field, 2 )
% isequal(f1, f2)
% % sum over first fiber dimension 
% f1 = std( field,1 )
%--------------------------------------------------------------------------
% Author: Fabian Telschow
%--------------------------------------------------------------------------
if nargin == 1
   if isa( varargin{1}, "Field")
       out = varargin{1};
       out.field = std( out.field, 0, out.D + out.fiberD );
   end
elseif nargin == 2
   if isa( varargin{1}, "Field")
       D = varargin{2};
       if isnumeric( D )
          if D > 0
            out = varargin{1};
            out.field = std( out.field, 0, out.D + D );
          else
            error( "Second inout must be a positive integer." )
          end
       end
   else
       error( "First input must be an object of class Field." )
   end
elseif nargin == 3
   if isa( varargin{1}, "Field")
       D = varargin{2};
       if isnumeric( D )
          if D > 0
            out = varargin{1};
            out.field = std( out.field, varargin{3}, out.D + D );
          else
            error( "Second inout must be a positive integer." )
          end
       end
   else
       error( "First input must be an object of class Field." )
   end
else
  error( "Up to two inputs are supported." )
end

end 