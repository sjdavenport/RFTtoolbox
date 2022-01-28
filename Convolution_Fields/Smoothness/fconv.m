function [ smoothed_data, ss ] = fconv( data, sep_kern, D, truncation, dx,...
                                         adjust_kernel )
% FCONV( data, sep_kern, D, truncation, dx, adjust_kernel ) provides a fast
% implementation for smoothing data using a separable kernel (e.g. an
% isotropic Gaussian kernel).
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  data        a Dim by nsubj array of data
%  sep_kern    possible options are
%               - a function handle giving a univariate kernel. This kernel
%                 is used for smoothing in all dimensions.
%               - a 1 by D cell array containing in each component the
%                 function handle for the univariate kernel to smooth the
%                 d-th dimension.
%               - a numeric or 1 by D vector, in which case fconv smoothes
%                 with a multivariate Gaussian kernel of FWHM = sep_kernel
%               - a SepKernel object, in which case the 'kernel' field is
%                 used for smoothing.
% Optional
%  D           the dimension of the input data, if not specified then the
%              data is smoothed assuming that there is only one realisation
%  truncation  either a numeric, a 1 x D vector or a 2 x D array containing
%              the value for the truncation of the kernel in each dimension.
%              If it is a vector the truncation is assumed to be symmetric
%              around the origin, i.e., in the dth dimension the vector the
%              kernel is evaluated on is -truncation(d):dx(d):truncation(d).
%              If it is a numeric the kernel in each dimension is evaluated
%              on the vector -truncation:dx:truncation.
%  dx          either a numeric or a vector which defines the stepsize for
%              the evaluation grid of the kernel, see also truncation
%              description. Default is 1.
%  adjust_kernel  1 by D numeric vector computing the smoothing kernel at
%                 an offset, i.e. K( x + adjust_kernel ) instead of K( x ).
%--------------------------------------------------------------------------
% OUTPUT
% smoothed_data     the smoothed data
% ss                the sum of squares of the kernel (useful for ensuring
%                   variance 1 isotropic fields
%--------------------------------------------------------------------------
% DEVELOPER TODOS
% Need to allow calculation of    ss (the sum of squares of the kernel for use
% in noisegen.m
%--------------------------------------------------------------------------
% EXAMPLES
% %% %% 1D Examples
% %% Simple 1D example
% lat_data = normrnd(0,1,1,100); FWHM = 3;
% smoothed_fconv = fconv(lat_data, FWHM);
% smoothed_spm = spm_conv(lat_data,FWHM);
% figure, clf,
% plot(smoothed_spm); hold on; plot(smoothed_fconv, '--')
% legend('spm\_conv', 'fconv') 
% %% Using SepKernel object
% sepK = SepKernel(1, 3);
% lat_data = normrnd(0,1,1,100);
% smoothed_fconv = fconv(lat_data, sepK);
% %% 1D multiple subjects
% nvox = 100; nsubj = 2; FWHM = 3; D = 1;
% lat_data = normrnd( 0, 1, nvox, nsubj );
% smoothed_fconv = fconv( lat_data, FWHM, D );
% smoothed_spm = zeros( nvox, nsubj );
% for n = 1:nsubj
%     smoothed_spm(:,n) = spm_conv( lat_data(:,n), FWHM );
% end
% plot( smoothed_spm, 'color',  [ 0 0.447 0.7410 ] ); hold on;
% plot( smoothed_fconv, '--', 'color', [ 0.85 0.325 0.0980 ] );
% 
% %% %% 2D Examples
% %% % Simple 2D example using numeric FWHM input
% lat_data = normrnd( 0, 1, 25, 25 ); FWHM = 5;
% smoothed_fconv = fconv( lat_data, FWHM ); 
% smoothed_spm = spm_conv( lat_data, FWHM );
% subplot(2,1,1)
% surf(smoothed_fconv)
% title('fconv')
% subplot(2,1,2)
% surf(smoothed_fconv)
% title('SPM\_conv')
% 
% %% % Simple 2D example using numeric vector FWHM input
% % Generate lattice data
% lat_data = normrnd( 0, 1, 25, 25 );
% % FWHM  for smoothing in each direction
% FWHM = [3, 7];
% % change the frequency where the kernel is evaluated. Voxels are considered
% % dx units away
% dx = 0.5;
% % Smooth the data
% smoothed_fconv = fconv( lat_data, FWHM ); 
% figure, clf,
% surf(smoothed_fconv)
% title('fconv defaults')
% 
% %% % Same result but filling the arguments manually
% % Create a Sep_Kernel objcet representing the Gaussian kernel
% K = SepKernel( 2, FWHM )
% % Smooth the data
% smoothed_fconv = fconv( lat_data, K.kernel, 2, K.truncation, 1, K.adjust );
% figure, clf,
% surf(smoothed_fconv)
% title('fconv manually filled')
% 
% %% % Change of size of voxels
% % change the frequency where the kernel is evaluated. Voxels are considered
% % dx units away
% dx = 0.5;
% % Smooth the data
% smoothed_fconv = fconv( lat_data, K.kernel, 2, K.truncation, dx, K.adjust );
% figure, clf,
% surf(smoothed_fconv)
% title('fconv voxels changed to one half')
% 
% %% %% 3D Examples
% Dim = [50,50,50]; lat_data = normrnd(0,1,Dim); halfDim = Dim(1)/2;
% FWHM = 3; D = 3;
% smoothed_spm = zeros(Dim);
% spm_smooth(lat_data, smoothed_spm, FWHM);
% smoothed_fconv = fconv(lat_data, FWHM);
% sigma = FWHM2sigma(FWHM); truncation = ceil(6*sigma);
% smoothed_fconv_spmkern = fconv(lat_data, @(x) spm_smoothkern(FWHM, x), D, truncation );
% smoothed_cfield = convfield( lat_data, FWHM, 0, D );
% figure, clf,
% plot(1:Dim(1),smoothed_fconv(:,halfDim,halfDim))
% hold on 
% plot(1:Dim(1),smoothed_spm(:,halfDim,halfDim))
% plot(1:Dim(1),smoothed_cfield(:,halfDim,halfDim), '--')
% plot(1:Dim(1),smoothed_fconv_spmkern(:,halfDim,halfDim), '--')
% legend('fconv', 'SPM', 'convfield', 'fconv_smoothkern')
% 
% figure, clf,
% plot(-truncation:truncation, spm_smoothkern(FWHM, -truncation:truncation))
% hold on
% plot(-truncation:truncation, GkerMV(-truncation:truncation, FWHM))
% title("Difference in Gaussian kernel and spm12 kernel")
% 
% % Compare speed to spm_smooth (much faster)
% Dim = [50,50,50]; lat_data = normrnd(0,1,Dim);
% tic; fconv(lat_data, FWHM); toc
% tic; smoothed_spm = zeros(Dim);
% tt = spm_smooth_mod(lat_data, smoothed_spm, FWHM); toc
% 
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport, Fabian Telschow
%--------------------------------------------------------------------------


%% Check mandatory input and get important constants
%--------------------------------------------------------------------------

% Find the dimension, of the data if not provided. Here it is assumed that
% nsubj = 1.
s_data = size( data );

% If D is not provided the data is smoothed assuming that nsubj = 1
if ~exist( 'D', 'var' )
    s_data = size( data );
    D = length( s_data );
    if D == 2
        if s_data(1) == 1 
            D = 1;
        elseif s_data(2) == 1
            D = 1;
            data = data';
        end
    end
end

% Set horizontal to vertical switch
if D == 1 && s_data(1) == 1
    horz2vert = 1;
else
    horz2vert = 0;
end

% Reject if D is to large.
if D > 3
    error( 'fconv not coded for dimension > 3' );
end

%%% Check sep_kernel input
% Make Kernel a cell with equal kernel functions
Kernel = cell( [ 1 D ] );

if isnumeric( sep_kern )
    % Numeric sep_kernel defines the FWHM
    sep_kern = SepKernel( D, sep_kern );

    for d = 1:D
        Kernel{d} = sep_kern.kernel{d};
    end
    
    % Get the standard truncation for the Gaussian kernel
    truncation = sep_kern.truncation;
    
elseif isa( sep_kern, 'SepKernel' )
    if length( sep_kern.kernel ) == D
        for d = 1:D
            Kernel{d} = sep_kern.kernel{d};
        end
        
        % Define the truncation and adjust_kernel, if not provided
        if ~exist( 'truncation', 'var' )
            if ~any( isnan( sep_kern.truncation ) )
                truncation = sep_kern.truncation;
            else
                error( "You need to provide a valid truncation." )
            end
        end
        
        if ~exist( 'adjust_kernel', 'var' )
            adjust_kernel = sep_kern.adjust;
        end
        
    else
        error( 'Your SepKernel object needs to have the same dimension as your data.' )
    end
else
    if iscell( sep_kern ) 
        if length( sep_kern ) == 1
            for d = 1:D
                Kernel{d} = sep_kern{1};
            end
        elseif length( sep_kern ) == D
            Kernel = sep_kern;
        else
            error( strcat( 'A Kernel must be specified for each ',...
                       'direction or 1 that will be the same for all' ) )
        end
    else
        for d = 1:D
            Kernel{d} = sep_kern;
        end
    end
    
    % If kernel is specified manually a truncation is mandatory
    if ~exist( 'truncation', 'var' )
        error( 'Need to specify truncation' )
    end
end

%% Check/add optional input
%--------------------------------------------------------------------------

%%% Check/transform truncation input
if numel( truncation(:) ) == 1  
    % Make truncation a 2 x D array.
    truncation = truncation * ones( [ 2 D ] );
    
elseif numel( truncation(:) ) == D
    % Make truncation a 2 x D array.
    if size( truncation, 1 ) == 1
        truncation = [ truncation; truncation ];
    else
        truncation = [ truncation, truncation ]';
    end
    
elseif all( size( truncation ) ~= [ 2 D ] )
    error( strcat( "If you want to specify truncation in",...
                   " every dimension then truncation must be 2 x D." ) )
end

%%% Check/transform dx input
if ~exist( 'dx', 'var' )
    % Set dx default.
    dx = ones( [ 1 D ] );
    
elseif numel( dx(:) ) == 1
    % Make dx a vector if numeric.
    dx = dx * ones( [ 1 D ] );
    
end

%%% Check adjust_kernel input
if ~exist( 'adjust_kernel', 'var' ) || isnan( sum( adjust_kernel(:) ) )
    % set default value
    adjust_kernel = zeros( 1, D );
    
elseif length( adjust_kernel ) ~= D
    error( strcat( 'The kernel adjustment must be of the ',...
                   'same dimension as the data, ensure that length(FWHM) = D!' ) )
end

%% Main function
%--------------------------------------------------------------------------

%%% If there are multiple subjects run fconv on each of them
if ( D < length( s_data ) && D > 1 ) ...
                                  || ( D == 1 && all( s_data > [ 1, 1 ] ) )

    % Preallocate the smoothed_data field
    smoothed_data = zeros( s_data );
    % Get the number of subjects
    index  = repmat( {':'}, 1, D );
    
    % Loop over subjects
    for J = 1:s_data( end )
        [tmp, ss] = fconv( squeeze( data( index{:}, J ) ),...
                                            Kernel, D, truncation, dx,...
                                            adjust_kernel );
        smoothed_data( index{:}, J ) = tmp;
    end
    
    return
end

%%% fconv for a single subject
% Preallocate a cell array for the truncation vectors
truncation_vecs = cell( [ 1 D ] );

for d = 1:D
    % Specify the values at which to evaluate the dth kernels at
    truncation_vecs{d} = ( -truncation( 1, d ):dx(d):truncation( 2, d ) )...
                         + adjust_kernel( d );
                     
end


% Calculate the kernel at the values in truncation vector
kernel_in_direction = cell( [ 1 D ] );
for d = 1:D
    % Get the kernel in the dth direction evaluate on truncation
    kernel_in_direction{d} = Kernel{d}( truncation_vecs{d} );
end


%%% Main loop, running convolution with a separable kernel in each direction
if D == 1
    if horz2vert
        data = data';
    end
    
    smoothed_data = conv( data, kernel_in_direction{1}', 'same' );
    
    if horz2vert
        smoothed_data = smoothed_data';
    end
    
    % Calculates the sum of the squares of the kernel
    ss = sum( kernel_in_direction{1}.^2 ) * dx;
    
elseif D == 2
    
    % The kernel in the x direction (x corresponds to the rows of the matrix)
    xside = kernel_in_direction{1}';
    % The kernel in the y direction (y corresponds to the rows of the matrix)
    yside = kernel_in_direction{2};
    
    % Smooth in the x direction
    smoothed_data = convn( data, xside, 'same' );
    % Smooth in the y direction
    smoothed_data = convn( smoothed_data, yside, 'same' );
                      
    % This is the sum of the squares of the kernel in the truncation box,
    % since the kernel is by assumption separable.
     ss = sum(sum( (kernel_in_direction{1}' * kernel_in_direction{2}).^2 )) * dx(1)*dx(2);
    
elseif D == 3
    
    % The kernel in the x direction
    xside = kernel_in_direction{1}';
    % The kernel in the y direction
    yside = kernel_in_direction{2};
    % The kernel in the z direction
    zside = zeros( 1, 1, length( kernel_in_direction{3} ) );
    zside( 1, 1, : ) = kernel_in_direction{3};
    
    % Smooth in the x direction
    smoothed_data = convn( data, xside, 'same' );
    % Smooth in the y direction
    smoothed_data = convn( smoothed_data, yside, 'same' );
    % Smooth in the z direction
    smoothed_data = convn( smoothed_data, zside, 'same' );
    
    % This is the sum of the squares of the kernel in the truncation box,
    % since the kernel is by assumption separable.
    [X,Y,Z] = meshgrid( kernel_in_direction{1}', kernel_in_direction{2}', ...
                        kernel_in_direction{3}' );
    ss = ( X .* Y .* Z ).^2;
    ss = sum(ss(:)) * dx(1)*dx(2)*dx(3);
    
else
    error('fconv not coded for dimension > 3')
end

end
