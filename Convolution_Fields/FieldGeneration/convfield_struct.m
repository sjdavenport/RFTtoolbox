function [ smooth_data, xvals_vecs ] = convfield_struct( lat_data, Kernel, resAdd,...
                        D, derivtype, enlarge)
% CONVFIELD( lat_data, Kernel, resAdd, D, derivtype ) generates a
% convolution field evaluated on an equidistant grid with resolution
% increased by adding resAdd voxels inbetween each voxel in each dimension.
% The field is derived from lattice data smoothed with an seperable kernel
% which can be specified.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  lat_data  a Dim by nsubj array of data
%  Kernel    either a structure or a numeric.
%            If structure it MUST contain the following fields:
%            - kernel:
%              if derivtype = 0: function handle computing the kernel.
%              if derivtype = 1: D by 1 cell array containing the function
%                                handles for the kernels of the derivative
%                                in each direction.
%              Optional if the field 'kernel' is numeric the convolution
%              field is assumed to be derived from a Gaussian kernel and 
%              the field 'truncation' if missing is chosen automatically.
%            - truncation:
%              a window around the points at which to evaluate the kernel
%              setting this parameter allows for quicker computation of the
%              convolution field and has very litle effect on the values of
%              the field for kernels that have light tails such as the 
%              Gaussian kernel.
%            The following fields are OPTIONAL:
%            - adjust_kernel:
%               a D by 1 vector that allows you to compute the convolution
%               field at an offset. This is useful for computing 
%               derivatives numerically.
%               Default is not to use this feature.
%
%            If Kernel is numeric, the convolution field is assumed to be
%            derived from the Gaussian kernel with FWHM = Kernel.
%            Truncation and adjust_kernel are set to be default values.
% Optional
% resAdd     the amount of voxels added equidistantly inbetween the
%            existing voxels. Default is 1.
% D          the dimension of the data, if this is left blank it is
%            assumed that nsubj = 1 and that the convolution field has 
%            the same numberof dimesions as lat_data
% derivtype  0/1/2, 0 returns the convolution field, 1 it's derivative
%            and 2 it's second derivative (at all points). Default is 0
%            i.e to return the field! Note if D > 1 then  derivtype = 2
%            has not been implemented yet.
% enlarge    a numeric which must be a positive integer or zero. The
%            convolution field is computed on a domain enlarged in each
%            direction by 'enlarge' voxels. Note if resAdd ~=0 the voxels
%            are in high resolution.
%            Default 0.
%--------------------------------------------------------------------------
% OUTPUT
% %%% 1D
% %% Smoothing with increased resolution
% nvox = 100; xvals = 1:nvox;
% xvals_fine = 1:0.01:nvox;
% FWHM = 3;
% lat_data = normrnd(0,1,1,nvox);
% cfield = inter_conv1D( lat_data, FWHM, 0.01);
% plot(xvals_fine,cfield)
% hold on
% smooth_data = convfield( lat_data, FWHM, 0.01, 1);
% plot(xvals_fine,smooth_data + 0.5)
% 
% %% Smoothing with the same resolution
% lat_data = normrnd(0,1,1,nvox);
% cfield = spm_conv(lat_data, FWHM);
% plot(xvals,cfield)
% hold on
% smooth_data = convfield( lat_data', FWHM, 1, 1 );
% plot(xvals, smooth_data)
% 
% %% Multiple subjects
% nsubj = 3; nvox = 100;
% lat_data = normrnd(0,1,nvox,nsubj);
% cfield = spm_conv(lat_data(:,1), FWHM);
% plot(1:nvox,cfield)
% hold on
% smooth_data = convfield( lat_data, FWHM, 1, 1 );
% plot(1:nvox,smooth_data(:,1))
% 
% %% 1D derivatives
% nvox = 100; h = 0.01; xvals = 1:nvox; xvals_fine = 1:0.01:nvox;
% lat_data = normrnd(0,1,nvox,1);
% smoothedfield = convfield( lat_data, FWHM, h, 1);
% deriv1 = convfield( lat_data, FWHM, h, 1, 1 );
% deriv2 = diff(smoothedfield)/h;
% plot(xvals_fine, deriv1 + 0.5)
% hold on 
% plot(xvals_fine(1:end-1), deriv2)
% 
% % 1D derivative (multiple subjects)
% nvox = 100; FWHM = 3; xvals = 1:nvox; lat_data = normrnd(0,1,1,nvox);
% aderiv = @(x) applyconvfield( x, lat_data, @(y) GkerMVderiv(y, FWHM) );
% deriv = convfield( lat_data, FWHM, h, 1, 1 );
% deriv(1), aderiv(1)
% 
% %% 2D
% Dim = [25,25];
% lat_data = normrnd(0,1,Dim);
% spmlatfield = spm_conv(lat_data, FWHM);
% smooth_data = convfield( lat_data, FWHM, 1, 2); %Same as spm_conv (except for boundary differences)
% fine_data = convfield( lat_data, FWHM, 0.25, 2); %Convolution eval
% 
% zlimits = [min(fine_data(:))-0.1, max(fine_data(:))+0.1];
% 
% subplot(1,2,1)
% surf(smooth_data)
% zlim(zlimits)
% title('Lattice Evaluation')
% subplot(1,2,2)
% surf(fine_data)
% zlim(zlimits)
% title('Convolution Field')
% 
% %% Matching to applyconvfield
% cfield = @(x) applyconvfield(x, lat_data, FWHM);
% smooth_data(20,20)
% cfield([20,20]')
% smooth_data(1,10)
% cfield([1,10]')
% 
% %% 2D derivatives
% Dim = [25,25];
% lat_data = normrnd(0,1,Dim); spacing = 1;
% derivfield = convfield( lat_data, FWHM, spacing, 2, 1);
% surf(reshape(derivfield(1,:), spacep(Dim,spacing)))
% title('2D 1st partial derivative of the convolution field')
% 
% %% Showing that the derivatives work
% Dim = [5,5];
% lat_data = normrnd(0,1,Dim);
% 
% point = [3,3]';
% spacing = 0.01;
% 
% spaced_point = spacep(point,spacing);
% derivfield = convfield( lat_data, FWHM, spacing, 2, 1);
% 
% % convolution derivatives
% derivfield(:,spaced_point(1), spaced_point(2))
% 
% % Derivative using applyconvfield
% aderiv = @(x) applyconvfield( x, lat_data, @(y) GkerMVderiv(y, FWHM)  );
% aderiv(point)
% 
% % Illustration on a fine lattice (not to be used in practice)
% smoothfield100 = convfield_dep( lat_data, FWHM, spacing, 2, 0);
% partialderiv_finelat(1) = (smoothfield100(spaced_point(1)+1, spaced_point(2)) - smoothfield100(spaced_point(1),spaced_point(2)))/spacing;
% partialderiv_finelat(2) = (smoothfield100(spaced_point(1), spaced_point(2) + 1) - smoothfield100(spaced_point(1),spaced_point(2)))/spacing;
% partialderiv_finelat'
% % note that it doesn't match perfectly because it's still a discrete
% % approximation, but that's why we want to use derivfield in the first
% % place!
% 
% % SPM (i.e. lattice) estimates of the derivative (quite off!)
% smoothfield_spm = convfield( lat_data, FWHM, 1, 2, 0);
% spm_derivs(1) = (smoothfield_spm(point(1)+1, point(2)) - smoothfield_spm(point(1),point(2)));
% spm_derivs(2) = (smoothfield_spm(point(1), point(2) + 1) - smoothfield_spm(point(1),point(2)));
% spm_derivs'
% 
% %% 2D derivatives (multiple subjects)
% Dim = [5,5]; nsubj = 20;
% lat_data = normrnd(0,1,[Dim, nsubj]);
% derivfield = convfield( lat_data, FWHM, 1, 2, 1)
% 
% %% 3D
% %Compare to SPM
% Dim = [10,10,10];
% FWHM = 3;
% lat_data = normrnd(0,1,Dim);
% subplot(1,2,1)
% spm_smooth_field = zeros(Dim); 
% spm_smooth(lat_data, spm_smooth_field, FWHM)
% surf(spm_smooth_field(:,:,5));
% title('Lattice Eval')
% subplot(1,2,2)
% cfield = convfield( lat_data, FWHM, 1, 3); %Convolution eval
% surf(cfield(:,:,5))
% title('Convolution Field Eval (no smoothing)')
% 
% %% Fine evaluation
% Dim = [10,10,10];
% spacing = 0.1; D = length(Dim); FWHM = 3; 
% slice = Dim(end)/2; fine_slice = Dim(end)/2/spacing;
% lat_data = normrnd(0,1,Dim);
% subplot(1,3,1)
% spm_smooth_field = zeros(Dim); 
% spm_smooth(lat_data, spm_smooth_field, FWHM)
% surf(spm_smooth_field(:,:,slice));
% title('Lattice Eval')
% subplot(1,3,2)
% cfield = convfield( lat_data, FWHM, spacing, D, 0, 0); %Convolution eval
% surf(cfield(:,:,fine_slice))
% title('Convolution Field Eval (Convn)')
% subplot(1,3,3)
% cfield_withspm = convfield( lat_data, FWHM, spacing, D, 0, 1); %Convolution eval (using spm_smooth)
% surf(cfield_withspm(:,:,fine_slice))
% title('Convolution Field Eval (SPM\_smooth)')
% 
% 
% %% Compare to applyconvfield
% lat_data = normrnd(0,1,Dim);
% acfield = @(x) applyconvfield(x, lat_data, FWHM);
% Dim = [10,10,10];
% spacing = 1; D = length(Dim); FWHM = 3; 
% cfield = convfield( lat_data, FWHM, 1, D, 0, 0); 
% cfield_withspm = convfield( lat_data, FWHM, 1, D, 0, 1);
% acfield([5,5,5]')
% cfield(5,5,5)
% cfield_withspm(5,5,5) %Even within the image the spm_smooth version is
% % off ( as SPM does something weird in the z-direction)
% acfield([1,1,10]')
% cfield(1,1,10)
% cfield_withspm(1,1,10) %SPM_smooth version is quite off on the boundary!
% 
% % Note that the differneces in spm_smooth appear to go away as the FWHM or the
% % spacing increases, so the spm_smooth version does well to evaluate ifne
% % convolution fields and does so very efficiently. So the difference could 
% % be caused by a truncation issue in spm_smooth?
% 
% %% % 3D derivatives (1 subject)
% Dim = [5,5,5]; D = length(Dim); FWHM = 3;
% lat_data = normrnd(0,1,Dim);
% derivfield = convfield( lat_data, FWHM, 1, D, 1)
% aderiv = @(x) applyconvfield( x, lat_data, @(y) GkerMVderiv(y, FWHM)  )
% dfeval = derivfield(:,3,3,3)
% aceval = aderiv([3,3,3]')
% spacing = 0.05;
% spaced_point = spacep( [3,3,3]', spacing);
% cfield_fine = convfield( lat_data, FWHM, spacing, D);
% pointeval = cfield_fine(spaced_point(1),spaced_point(2),spaced_point(3));
% plusxeval = cfield_fine(spaced_point(1)+1,spaced_point(2),spaced_point(3));
% plusyeval = cfield_fine(spaced_point(1),spaced_point(2)+1,spaced_point(3));
% pluszeval = cfield_fine(spaced_point(1),spaced_point(2),spaced_point(3)+1);
% derivx = (plusxeval - pointeval)/spacing
% derivy = (plusyeval - pointeval)/spacing
% derivz = (pluszeval - pointeval)/spacing
% 
% %% 3D derivatives (Multiple subjects)
% Dim = [5,5,5]; D = length(Dim); FWHM = 3; nsubj = 2;
% lat_data = normrnd(0,1,[Dim,nsubj]);
% derivfields = convfield( lat_data, FWHM, 1, D, 1)
% aderiv = @(x) applyconvfield( x, lat_data(:,:,:,2), @(y) GkerMVderiv(y, FWHM)  )
% dfeval = derivfields(:,3,3,3,2)
% aceval = aderiv([3,3,3]')
% 
% %% Adjusting the field (1D)
% nvox = 10; D = 1; FWHM = 2; lat_data = normrnd(0,1,1,nvox);
% spacing = 0.1; adjust_kernel = 0.1';
% [smoothfield, xvals_vecs] =  convfield( lat_data, FWHM, spacing, D, 0, -1);
% 
% [adjust_field, xvals_vecs_adjust] = convfield( lat_data, FWHM, 1, D, 0, -1, adjust_kernel );
% 
% plot(xvals_vecs{1}, smoothfield)
% hold on
% plot(xvals_vecs_adjust{1},adjust_field)
% 
% acfield = @(tval) applyconvfield(tval, lat_data, FWHM);
% adjust_field(3)
% acfield(3.1)
% 
% %% Adjusting the field (3D)
% Dim = [10,10,10]; D = length(Dim); FWHM = 1.5;
% lat_data = normrnd(0,1,Dim); spacing = 0.1;
% [smoothfield, xvals_vecs] = convfield( lat_data, FWHM, spacing, D, 0, -1);
% 
% adjust_kernel = [0.1,0,0]';
% [adjust_field, xvals_vecs_adjust] = convfield( lat_data, FWHM, 1, D, 0, -1, adjust_kernel );
% 
% point = [1.1,1,1]'; spaced_point = spacep(point, spacing);
% plot(xvals_vecs{1}, smoothfield(:,spaced_point(2), spaced_point(3)))
% hold on
% plot(xvals_vecs_adjust{1},adjust_field(:,point(2),point(3)))
% acfield = @(tval) applyconvfield(tval, lat_data, FWHM) 
% adjust_field(1,1)
% acfield(point)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport, Fabian Telschow                                              
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Get important constants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% get constants from the input lat_data field
% get size of the input data
slatdata = size( lat_data );
% get number of dimensions of input data
D_latdata = length( slatdata );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% add/check optional values and define the kernel structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% D input
% If no dimension is specfied it is assumed that nsubj = 1 and you just
% want to smooth a single field
if ~exist( 'D', 'var' )
    D = D_latdata;
end

% get size of the domain and number of subjects from lat_data
if D > 1
    if D_latdata == D
        nsubj = 1;
        Dim   = slatdata;
    else
        nsubj = slatdata( end );
        Dim   = slatdata( 1 : end-1 );
    end
else
    % This loop is included to find Dim in 1D and to allow for different 
    % types of 1D data input
    % A parameter that will determine whether to transpose the data or not
    vert2horz = 0;
    if D_latdata == 2 && slatdata(1) == 1
        nsubj     = 1;
        Dim       = slatdata( slatdata > 1 );
        vert2horz = 1; % I.e. if a horizontal field is entered it returns one as well
    elseif D_latdata == 2 && slatdata(2) == 1
        nsubj = 1;
        Dim   = slatdata( slatdata > 1 );
    else %I.e. if D_lat_data == 2 and slatdata(1) and slatdata(2) are both > 1
        % or if D_lat_data > 2 and there are some 1 dimensional dimensions
        % hanging around for some reason.
        slatdata_squeezed = size( squeeze( lat_data ) );
        Dim   = slatdata_squeezed(1);
        nsubj = slatdata_squeezed(2);
    end
end

% index to easier code multidimensional
indexD = repmat( {':'}, 1, D );

%%%%%%% resAdd input
% add resolution increase if missing
if ~exist( 'resAdd', 'var' )
    resAdd = 1;
end
% get the difference between voxels with resolution increase
dx = 1 / ( resAdd + 1 );

% !!!small fix, such that old code still works, need to be removed at some
% point
if resAdd < 1 && 0 < resAdd
    resAdd = 1;
end

% reject input, if resAdd is to large in 3D
if D == 3 && ( resAdd > 18 )
    error( 'In 3D you shouldn''t use such high resolution for memory reasons' )
end

%%%%%%% derivtype input
% add derivative type if missing
if ~exist( 'derivtype', 'var' )
    derivtype = 0;
end


%%%%%%% enlarge input
if ~exist( 'enlarge', 'var' )
    enlarge = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% check and prepare the Kernel input structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% indicating whether kernel needs to be adjusted. Will be changed to 1 if
% needed
use_adjust = 0;

if isnumeric( Kernel ) % If Kernel is numeric use Gaussian Kernel
    % set the FWHM parameter to be the user input adjusted by the
    % resolution factor
    FWHM = Kernel  / dx;
    % change numerical Kernel input to structure
    Kernel = struct();
    % Obtain the parameter of kernel in form of the FWHM adjusted for
    % the dx
    sigma = FWHM2sigma( FWHM );
    % Set default truncation
    Kernel.truncation = ceil( 4*sigma );
    % Set default adjust_kernel
    Kernel.adjust_kernel = zeros( D, 1 );
    
    % Default kernel is isotropic Gaussian kernel, respectively its
    % derivatives
    if derivtype == 0
        Kernel.kernel = @(x) Gker( x, FWHM );
    elseif derivtype == 1
        switch D
            case 1
                Kernel.kernel =  { @(x) Gkerderiv( x, FWHM ) };
            case 2
                Kernel.kernel = cell([ 1 D ]);
                Kernel.kernel{1,1} = { @(x) Gkerderiv( x, FWHM ),...
                                       @(y) Gker( y, FWHM ) };
                Kernel.kernel{1,2} = { @(x) Gker( x, FWHM ),...
                                       @(y) Gkerderiv( y, FWHM ) };
            case 3
                Kernel.kernel = cell([ 1 D ]);
                Kernel.kernel{1,1} = { @(x) Gkerderiv( x, FWHM ),...
                                       @(y) Gker( y, FWHM ),...
                                       @(z) Gker( z, FWHM ) };
                Kernel.kernel{1,2} = { @(x) Gker( x, FWHM ),...
                                       @(y) Gkerderiv( y, FWHM ),...
                                       @(z) Gker( z, FWHM ) };
                Kernel.kernel{1,3} = { @(x) Gker( x, FWHM ),...
                                       @(y) Gker( y, FWHM ),...
                                       @(z) Gkerderiv( z, FWHM ) };
        end
    else
        error( 'This setting has not been coded yet' )
    end
    
elseif isstruct( Kernel )
    if ~isfield( Kernel, 'kernel' )
        error( strcat( "The structure 'Kernel' needs to have at",...
                       " least a field called 'kernel'!\n%s" ),...
               "Please read the input instructions!" )
    else
        %%%% set default of field 'adjust_kernel' if missing 
        if ~isfield( Kernel, 'adjust_kernel' )
            Kernel.adjust_kernel = zeros( D, 1 );            
        else
            % Check adjust kernel field
            if any( Kernel.adjust_kernel )
                use_adjust = 1;
            end
            % make sure adjust_kernel hat the correct orientation
            if size( Kernel.adjust_kernel, 1 ) == 1
                Kernel.adjust_kernel = Kernel.adjust_kernel';
            end
            % ensure that adjust_kernel has D entries
            if length( Kernel.adjust_kernel ) ~= D
                error( 'The kernel adjustment must be of the right dimension' )
            end
        end
        %%%% check whether field 'kernel' is numeric or not
        if ~isnumeric( Kernel.kernel )
            % check whether the mandatory truncation field exists
            if ~isfield( Kernel, 'truncation' )
                error( strcat( "For a general Kernel you need to",...
                               " provide a field 'truncation'!\n%s" ),...
                       "Please read the input instructions!" )
            end
            
            % reject input if derivtype not equal to zero, since the
            % derivative cannot be computed in general and must be provided
            % by the user for a general kernel in form of the derivative of
            % the Kernel
            if derivtype ~= 0
                error( strcat( "To obtain the derivative of a convolution",...
                               " field from a general kernel you need to",...
                               " provide the derivative Kernel.\n%s" ),...
                       strcat( "Then, if you set derivtype = 0, this",...
                               " function computes your derivative!" ) )
            end
        else
            % set the FWHM parameter to be the user input adjusted by the
            % resolution factor
            FWHM = Kernel.kernel / dx;
            % Obtain the parameter of kernel in form of the FWHM adjusted
            % for the dx
            sigma = FWHM2sigma( FWHM );
            % Set default truncation if missing
            if ~isfield( Kernel, "truncation" )
                Kernel.truncation = ceil( 4*sigma );
            end
            % Set default adjust_kernel
            if ~isfield( Kernel, "adjust_kernel" )
                Kernel.adjust_kernel = zeros( D, 1 );
            end
            
            % Check adjust kernel field
            if any( Kernel.adjust_kernel )
                use_adjust = 1;
            end
            % make sure adjust_kernel hat the correct orientation
            if size( Kernel.adjust_kernel, 1 ) == 1
                Kernel.adjust_kernel = Kernel.adjust_kernel';
            end
            % ensure that adjust_kernel has D entries
            if length( Kernel.adjust_kernel ) ~= D
                error( 'The kernel adjustment must be of the right dimension' )
            end
            
            % Default kernel is isotropic Gaussian kernel, respectively its
            % derivatives
            if derivtype == 0
                Kernel.kernel = @(x) Gker( x, FWHM );
            elseif derivtype == 1
                switch D
                    case 1
                        Kernel.kernel = { @(x) Gkerderiv( x, FWHM ) };
                    case 2
                        Kernel.kernel = cell([ 1 D ]);
                        Kernel.kernel{1,1} = { @(x) Gkerderiv( x, FWHM ),...
                                               @(y) Gker( y, FWHM ) };
                        Kernel.kernel{1,2} = { @(x) Gker( x, FWHM ),...
                                               @(y) Gkerderiv( y, FWHM ) };
                    case 3
                        Kernel.kernel = cell([ 1 D ]);
                        Kernel.kernel{1,1} = { @(x) Gkerderiv( x, FWHM ),...
                                               @(y) Gker( y, FWHM ),...
                                               @(z) Gker( z, FWHM ) };
                        Kernel.kernel{1,2} = { @(x) Gker( x, FWHM ),...
                                               @(y) Gkerderiv( y, FWHM ),...
                                               @(z) Gker( z, FWHM ) };
                        Kernel.kernel{1,3} = { @(x) Gker( x, FWHM ),...
                                               @(y) Gker( y, FWHM ),...
                                               @(z) Gkerderiv( z, FWHM ) };
                end
            else
                error( 'This setting has not been coded yet' )
            end            
        end 
    end
end

%%%% move the kernel parameters into simpler variables such that not always
%%%% the structure needs to be used
adjust_kernel = Kernel.adjust_kernel;
truncation = Kernel.truncation;
Kernel = Kernel.kernel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% main part of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dimensions for domain of the field with increased resolution
Dimhr = ( Dim - 1 ) * resAdd + Dim; %Dimhr = Dim with high resolution
% modify Dimhr by adding enlarge voxels to all sides
if enlarge ~= 0 
    if D == 1
        Dimhr = Dimhr + [ 2 * enlarge, 0 ];
    else
        Dimhr = Dimhr + 2 * enlarge;
    end
end

% Setting up the default xvals_vecs
xvals_vecs  = cell( 1, D );
if use_adjust
    for d = 1:D
        xvals_vecs{d} = ( ( 1 - enlarge*dx ):dx:( Dim(d) + enlarge*dx ) )...
                            + adjust_kernel(d);
    end
else
    for d = 1:D
        xvals_vecs{d} = ( 1 - enlarge*dx ):dx:( Dim(d) + enlarge*dx );
    end
end

% create index to fill the original data at the correct voxels of the high
% resolution data array
index = cell( [ 1 D ] );
for d = 1:D
    index{d} = ( enlarge + 1 ):( resAdd + 1 ):( Dimhr(d) - enlarge );
end

% Increase the resolution of the raw data by introducing zeros
expanded_lat_data = zeros( [ Dimhr, nsubj] );
expanded_lat_data( index{:}, : ) = lat_data;

%%%% Main loop: Calculation of convolution fields
if D == 100 %@Sam: change back to 1 to see the error.
            %      I would very much like to remove the whole bit and always
            %      force the user to obey our convention of columns as
            %      samples!
    % Points at which to evaluate the Kernel
    gridside  = -truncation:dx:truncation;
     
    % Field adjustment if that is specified (note default is to bypass this loop)
    if use_adjust
        gridside = gridside + adjust_kernel;
%       gridside = fliplr(gridside - adjust_kernel);
    end
    
    if derivtype == 1
        Kernel = Kernel{1};
    end
    
    % Evaluates the kernel to get the filter for the convolution
    h = Kernel( gridside );
    % Performs convolution to get the convolution fields
    smooth_data = convn( expanded_lat_data', h, 'same' )';
    % @Sam: can you explain the factor
    smooth_data = smooth_data / dx^( D + derivtype );
    
    % Transpose the data to return a horizontal output, if initial input
    % was horizontal
    if vert2horz && nsubj == 1
        smooth_data = smooth_data';
    end
    
elseif D < 4
    % Run the smoothing using fconv
    if derivtype == 0
        % Calculates the convolution field
        smooth_data = fconv( expanded_lat_data, Kernel, D,...
                             truncation, adjust_kernel );
        % @Sam: can you explain the factor, naively I would have thought
        % that should be accounted for in the kernel and its derivative
        smooth_data = smooth_data / dx^D;
        
    elseif derivtype == 1
        % prelocate the output field for speed
        smooth_data = ones( [ Dimhr nsubj D ] );
        
        % Calculates the derivatives of the convolution field
        for d = 1:D
            smooth_data(indexD{:},:,d) = fconv( expanded_lat_data,...
                                                Kernel{d}, D, ...
                                                truncation );
        end
        % @Sam: can you explain the factor
        smooth_data = squeeze( smooth_data ) / dx^(D+1);
        
    else
        error( 'Higher derivatives are not supported' )  
    end
    
else
    error('D != 1,2,3 has not been implemented yet!')
end

return