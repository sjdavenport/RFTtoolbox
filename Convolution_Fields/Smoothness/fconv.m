function [ smoothed_data, ss ] = fconv( data, sep_kern, D, truncation, dx,...
                                         adjust_kernel )
% FCONV2( data, sep_kern, truncation, adjust_kernel ) provides a fast
% implementation for smoothing data using a separable kernel (e.g. an
% isotropic Gaussian kernel).
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  data        a Dim by nsubj array of data
%  sep_kern    a function handle giving a separable kernel. If this is 
%              instead numeric fconv smoothes with an isotropic Gaussian 
%              kernel with sep_kern as the FWHM (see EXAMPLES section)
% Optional
%  D           the dimension
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
%  adjust_kernel  @Sam please write this!
%--------------------------------------------------------------------------
% OUTPUT
% smoothed_data     the smoothed data
% ss                the sum of squares of the kernel (useful for ensuring
%                   variance 1 isotropic fields
%--------------------------------------------------------------------------
% EXAMPLES
% %1D
% lat_data = normrnd(0,1,1,100); FWHM = 3;
% smoothed_fconv = fconv(lat_data, FWHM);
% smoothed_spm = spm_conv(lat_data,FWHM);
% plot(smoothed_spm); hold on; plot(smoothed_fconv)
% legend('spm\_conv', 'fconv') 
% 
% %1D multiple subjects
% nvox = 100; nsubj = 2; lat_data = normrnd(0,1,nvox,nsubj); FWHM = 3; D = 1;
% smoothed_fconv = fconv(lat_data, FWHM, D)
% smoothed_spm = zeros(nvox, nsubj);
% for n = 1:nsubj
%     smoothed_spm(:,n) = spm_conv(lat_data(:,n),FWHM);
% end
% plot(smoothed_spm, 'color',[0.85 0.325 0.0980]); hold on;
% plot(smoothed_fconv, 'color',[0 0.447 0.7410]);
% 
% % 2D 
% lat_data = normrnd(0,1,25,25); FWHM = 3;
% smoothed_fconv = fconv(lat_data, FWHM); 
% smoothed_spm = spm_conv(lat_data,FWHM)
% subplot(2,1,1)
% surf(smoothed_fconv)
% title('fconv')
% subplot(2,1,2)
% surf(smoothed_fconv)
% title('SPM\_conv')
% 
% % 3D
% Dim = [50,50,50]; lat_data = normrnd(0,1,Dim); halfDim = Dim(1)/2;
% D = length(Dim); FWHM = 3; D = 3;
% smoothed_spm = zeros(Dim);
% spm_smooth(lat_data, smoothed_spm, FWHM);
% smoothed_fconv = fconv(lat_data, FWHM);
% sigma = FWHM2sigma(FWHM); truncation = ceil(6*sigma);
% smoothed_fconv_spmkern = fconv(lat_data, @(x) spm_smoothkern(FWHM, x), D, truncation );
% smoothed_cfield = convfield( lat_data, FWHM, 1, D, 0, 0);
% plot(1:Dim(1),smoothed_fconv(:,halfDim,halfDim))
% hold on 
% plot(1:Dim(1),smoothed_spm(:,halfDim,halfDim))
% plot(1:Dim(1),smoothed_cfield(:,halfDim,halfDim), '--')
% plot(1:Dim(1),smoothed_fconv_spmkern(:,halfDim,halfDim), '--')
% legend('fconv', 'SPM', 'convfield', 'fconv_smoothkern')
% 
% plot(-truncation:truncation, spm_smoothkern(FWHM, -truncation:truncation))
% hold on
% plot(-truncation:truncation, GkerMV(-truncation:truncation, FWHM))
% 
% % Compare speed to spm_smooth (much faster)
% Dim = [50,50,50]; lat_data = normrnd(0,1,Dim);
% tic; fconv(lat_data, FWHM); toc
% tic; smoothed_spm = zeros(Dim);
% tt = spm_smooth_mod(lat_data, smoothed_spm, FWHM); toc
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport, Fabian Telschow
%--------------------------------------------------------------------------


%% Check mandatory input and get important constants
%--------------------------------------------------------------------------

% Find the dimension, of the data if not provided. Here it is assumed that
% nsubj = 1.
s_data = size( data );
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

% Reject if D is to large.
if D > 3
    error( 'fconv not coded for dimension > 3' );
end

%%% sep_kernel input
% Make Kernel a cell with equal kernel functions
Kernel = cell( [ 1 D ] );

if isnumeric( sep_kern )
    % Numeric sep_kernel defines the FWHM
    FWHM   = sep_kern;

    for d = 1:D
        Kernel{d} = @(x) Gker( x, FWHM );
    end
    % Get the standard truncation for the Gaussian kernel
    sigma  = FWHM2sigma( FWHM );
    truncation = ceil( 4 * sigma );
    
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

%%% truncation input
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

%%% dx input
if ~exist( 'dx', 'var' )
    % Set dx default.
    dx = ones( [ 1 D ] );
    
elseif numel( dx(:) ) == 1
    % Make dx a vector if numeric.
    dx = dx * ones( [ 1 D ] );
    
end

%%% adjust_kernel input
if ~exist( 'adjust_kernel', 'var' ) || isnan( sum( adjust_kernel(:) ) )
    % set default value
    adjust_kernel = zeros( 1, D );
    
elseif length( adjust_kernel ) ~= D
    error( strcat( 'The kernel adjustment must be of the ',...
                   'same dimension as the data' ) )
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
        smoothed_data( index{:}, J ) = fconv( squeeze( data( index{:}, J ) ),...
                                            Kernel, D, truncation, dx,...
                                            adjust_kernel );
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
kernel_in_direction = cell( D );
for d = 1:D
    % Get the kernel in the dth direction evaluate on truncation
    kernel_in_direction{d} = Kernel{d}( truncation_vecs{d} );
end


%%% Main loop, running convolution with a separable kernel in each direction
if D == 1
    
    smoothed_data = conv( data, kernel_in_direction{1}', 'same' );
    
    % Calculates the sum of the squares of the kernel
    ss = sum( kernel_in_direction{1}.^2 );
    
elseif D == 2
    
    % The kernel in the x direction (x corresponds to the rows of the matrix)
    xside = kernel_in_direction{1}';
    % The kernel in the y direction (y corresponds to the rows of the matrix)
    yside = kernel_in_direction{2};
    
    % Smooth in the x direction
    smoothed_data = convn( data, xside, 'same' );
    % Smooth in the y direction
    smoothed_data = convn( smoothed_data, yside, 'same' );
%     smoothed_data = convn(udside, smoothed_data, 'same');  %Smooth in the up down direction
                      
    % This is the sum of the squares of the kernel in the truncation box,
    % since the kernel is by assumption separable.
%     ss = sum(sum( kernel_in_direction{1}' * kernel_in_direction{2} ));
    
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
%     tmp = kernel_in_direction{1}' * kernel_in_direction{2}
%     ss =  sum(sum( kernel_in_direction{1}' * kernel_in_direction{2} ));
    
else
    error('fconv not coded for dimension > 3')
end

end

% Deprecated
% if D == 1
%     smoothed_data = conv(smoothed_data, kernel_in_direction(d,:),'same'); 
%     % Could use convn here, but this way gives the same result no matter 
%     % whether data is a row or a column vector
% else
%     for d = 1:D
%         smoothed_data = convn(smoothed_data, kernel_in_direction(d,:),'same');
%     end
% end
