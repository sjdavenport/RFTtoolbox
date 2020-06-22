%% %% D = 2 
%% % General example
% generate a separable Gaussian kernel
gK = SepKernel( 2, [ 3, 6 ] );
% Define xvals vector
xvals = cell( [1 2] );
xvals{1} = -23:0.5:15;
xvals{2} = -5:0.25:3;
% evaluate the kernel and plot
evals = eval( gK, xvals, 0 );
figure, clf,
imagesc(evals)
%% % Default value for dx numeric
% evaluate the kernel and plot
evals = eval( gK, 1, 0 );
figure, clf,
imagesc( evals )
title( "dx=1" )
% evaluate the kernel and plot
figure, clf,
evals = eval( gK, 0.5, 0 );
imagesc( evals )
title( "dx=0.5" )
% evaluate the kernel and plot
figure, clf,
evals = eval( gK, 0.25, 0 );
imagesc( evals )
title( "dx=0.25" )