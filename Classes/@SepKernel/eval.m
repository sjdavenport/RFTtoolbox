function vals = eval( obj, xvals, derivtype )
% eval( obj ) evaluates a SepKernel object or its derivatives on
% a given grid.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  obj   an SepKernel object
%  xvals an 1 by D cell array containing on which grid values to
%        evaluate the d-th kernel
% Optional
%  derivtype possible values 0/1/2, 0 evaluates the seperable kernel, 1 its
%            derivatives and 2 its second derivatives. Default is 0.
%--------------------------------------------------------------------------
% OUTPUT
% vals  a obj.D by obj.D cell array containing an SepKernel object
%       for each element of the Hessian matrix of the inputed
%       SepKernel.
%
%--------------------------------------------------------------------------
% EXAMPLES
% %% D = 2 
% %% % General example
% % generate a separable Gaussian kernel
% gK = SepKernel( 2, [ 3, 6 ] );
% % Define xvals vector
% xvals = cell( [1 2] );
% xvals{1} = -23:0.5:15;
% xvals{2} = -5:0.25:3;
% % evaluate the kernel and plot
% evals = eval( gK, xvals, 0 );
% figure, clf,
% imagesc(evals)
% % evaluate the derivatives of the kernel and plot
% % Note that a 1 x 2 cell array is outputted containing
% % the evaluated kernel
% evals = eval( gK, xvals, 1 );
% figure, clf,
% imagesc(evals{1})
% title('x-direction derivative')
% figure, clf,
% imagesc(evals{2})
% title('y-direction derivative')
% %% % Default value for dx numeric
% % evaluate the kernel and plot
% evals = eval( gK, 1, 2 );
% figure, clf,
% imagesc( evals )
% title( "dx=1" )
% % evaluate the kernel and plot
% figure, clf,
% evals = eval( gK, 0.5, 0 );
% imagesc( evals )
% title( "dx=0.5" )
% % evaluate the kernel and plot
% figure, clf,
% evals = eval( gK, 0.25, 0 );
% imagesc( evals )
% title( "dx=0.25" )
%
%--------------------------------------------------------------------------
% Author: Fabian Telschow
%--------------------------------------------------------------------------

%% Check mandatory input
%--------------------------------------------------------------------------

if iscell( xvals )
  if length( xvals ) == 1
      tmp = xvals{1};
      xvals = cell( [ 1 obj.D ] );
      for d = 1:obj.D
          xvals{d} = tmp;
      end
  elseif length( xvals ) ~= obj.D
      error( "Please, provide xvalues for each dimension." )
  end
elseif numel( xvals ) == 1
  dx = xvals;
  xvals = cell( [ 1 obj.D ] );
  for d = 1:obj.D
      if all( size( obj.truncation ) == [ 2 obj.D ] )
          xvals{d} = -obj.truncation(1,d):dx:obj.truncation(2,d);
      elseif all( size( obj.truncation ) == [ 1 obj.D ] )
          xvals{d} = -obj.truncation(d):dx:obj.truncation(d);
      else
          error("Your truncation field in your SepKernel is incorrectly specified.")
      end
  end
else
  tmp = xvals;
  xvals = cell( [ 1 obj.D ] );
  for d = 1:obj.D
      xvals{d} = tmp;
  end
end

%% Main function
%--------------------------------------------------------------------------

switch obj.D
  case 1
      locs1 = xvals{1};
      locs{1} = locs1(:);
  case 2
      % Get the meshgrid
      [ locs1, locs2 ] = meshgrid( xvals{1}, xvals{2} );

      % Make the locations into a cell structure
      locs    = cell( [ 1 obj.D ] );
      locs{1} = locs1(:);
      locs{2} = locs2(:);

      % Get the size of the meshgrid
      slocs = size( locs1 );
  case 3
      % Get the meshgrid
      [ locs1, locs2, locs3 ] = meshgrid( xvals{1},...
                                          xvals{2}, xvals{3} );

      % Make the locations into a cell structure
      locs    = cell( [ 1 obj.D ] );
      locs{1} = locs1(:);
      locs{2} = locs2(:);
      locs{3} = locs3(:);

      % Get the size of the meshgrid
      slocs = size( locs1 );
end

% Evaluate the appropriate kernel on the grid
switch derivtype
  case 0
      % Evaluate the kernel on the meshgrid
      vals = ones( [ prod( slocs ) 1 ] );
      for d = 1:obj.D
          vals = vals .* obj.kernel{d}( locs{d}(:) );
      end

      % Shape output back toslocs output size
      if obj.D >1
        vals = reshape( vals, slocs );
      end
  case 1
      % Get the Gradient SepKernel objects of obj
      Dobj = Gradient( obj );
      % Initialize output cell
      vals = cell( [ 1 obj.D ] );
      
      for d = 1:obj.D
        vals{d} = eval( Dobj{d}, xvals, 0 );
      end
  case 2
      % Get the Hessian SepKernel objects of obj
      Dobj = Hessian( obj );
      % Initialize output cell
      vals = cell( [ obj.D obj.D ] );
      
      for d = 1:obj.D
          for dd = 1:obj.D
            vals{d, dd} = eval( Dobj{d, dd}, xvals, 0 );
          end
      end
end

return