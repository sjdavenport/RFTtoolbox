function [ smooth_data, xvals_vecs ] = convfield( lat_data, Kernel, resAdd,...
                        D, derivtype, truncation, adjust_kernel)
% CONVFIELD( lat_data, Kernel, resAdd, D, derivtype ) generates a
% convolution field evaluated on an equidistant grid with resolution
% increased by adding resAdd voxels inbetween each voxel in each dimension.
% The field is derived from lattice data smoothed with an seperable kernel
% which can be specified.
%--------------------------------------------------------------------------
% ARGUMENTS
% lat_data    a Dim by nsubj array of data
% FWHM        the FWHM of the kernel with which to do smoothing
% resAdd      the amount of voxels added equidistantly inbetween the
%             existing voxels. Default is 1.
% D           the dimension of the data, if this is left blank it is
%             assumed that nsubj = 1 and that the convolution field has 
%             the same numberof dimesions as lat_data
% derivtype   0/1/2, 0 returns the convolution field, 1 it's derivative
%             and 2 it's second derivative (at all points). Default is 0
%             i.e to return the field! Note if D > 1 then  derivtype = 2
%             has not been implemented yet.
% truncation  a window around the points at which to evaluate the kernel
%             setting this parameter allows for quicker computation of
%             the convolution field and has very litle effect on the values
%             of the field for kernels that have light tails such as the
%             Gaussian kernel. Default (which is recorded by setting
%             truncation = -1 or not including it) results in a truncation of 
%             4*FWHM2sigma(Kernel). 
% adjust_kernel    a D by 1 vector that allows you to compute the
%               convolution field at an offset. This is useful for computing 
%               derivatives numerically. Default is not to use this feature.
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
%%%% Check input and get important constants from the mandatory input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get size of the input data
slatdata = size( lat_data );
% get number of dimensions of input data
D_latdata = length( slatdata );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% add/check optional values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add spacing if missing
if nargin < 3
    resAdd = 1;
end

% Gives the difference between voxels with resolution increase
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

% add derivative type if missing
if nargin < 5
    derivtype = 0;
end

% If no dimension is specfied it is assumed that nsubj = 1 and you just
% want to smooth a single field
if nargin < 4
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

% Determine default kernel adjustments
if nargin < 7
    use_adjust = 0;
    adjust_kernel = zeros( D, 1 );
else
    if any( adjust_kernel )
        use_adjust = 1;
    end
    if size( adjust_kernel, 1 ) == 1
        adjust_kernel = adjust_kernel';
    end
    if length( adjust_kernel ) ~= D
        error( 'The kernel adjustment must be of the right dimension' )
    end
end

% Default Kernel
if isnumeric( Kernel )
    FWHM = Kernel;
    if nargin < 6
        truncation = -1;
    end

    if truncation == -1
        % Obtain the parameter of kernel in form of the FWHM adjusted for
        % the dx
        sigma = FWHM2sigma( FWHM / dx );
        % Set default truncation
        truncation = ceil( 4*sigma );
    end
    
    if derivtype == 0
        % Default kernel is isotropic Gaussian kernel
        Kernel = @(x) Gker( x, FWHM );
    elseif derivtype == 1
        Kernel =  @(x) Gkerderiv( x, FWHM );
    else
        error( 'This setting has not been coded yet' )
    end
else
    % Need a default for truncation for no gaussian kernels!
    if derivtype == 1
        error( 'derivtype > 0 for non-isotropic or non-Gaussian kernels has not been implemented yet' )
%         Kernel = @(x) getderivs( Kernel(x), D );
    elseif derivtype > 1
        error( 'derivtype > 1 for non-isotropic or non-Gaussian kernels has not been implemented yet' )
    end
end

% Setting up the default xvals_vecs
xvals_vecs  = cell( 1, D );
if use_adjust
    for d = 1:D
        xvals_vecs{d} = ( 1:dx:Dim(d) ) + adjust_kernel(d);
    end
else
    for d = 1:D
        xvals_vecs{d} = 1:dx:Dim(d);
    end
end

% Dimensions for domain of the field with increased resolution
Dimhr = ( Dim - 1 ) * resAdd + Dim; %Dimhr = Dim with high resolution

% Points at which to evaluate the Kernel
gridside  = -truncation:dx:truncation;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% main part of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main loop: Calculation of convolution fields
if D == 1
    % Increase the resolution of the raw data by introducing zeros
    expanded_lat_data = zeros( [ Dimhr, nsubj] );
    expanded_lat_data( 1:(resAdd + 1):end, : ) = lat_data;
     
    % Field adjustment if that is specified (note default is to bypass this loop)
    if use_adjust
        gridside = gridside + adjust_kernel;
%       gridside = fliplr(gridside - adjust_kernel);
    end
    
    h = Kernel(gridside); %Evaluates the kernel
    smooth_data = convn( expanded_lat_data', h, 'same' )'; % Performs convolution
    
    if vert2horz && nsubj == 1
        smooth_data = smooth_data'; % Transpose the data to return a horizontal output
    end
elseif D == 2
    % Increase the resolution of the raw data by introducing zeros
    expanded_lat_data = zeros( [ Dimhr, nsubj ] );
    expanded_lat_data( 1:( resAdd + 1 ):end, 1:( resAdd + 1 ):end, : ) = lat_data;
    
    % Adjust the FWHM to account for the dx
    adjusted_FWHM = FWHM / dx;
    
    % Run the smoothing
    if derivtype == 0 % Calculates the convolution field
        % Only set up for Gaussian kernels atm!
        smooth_data = fconv( expanded_lat_data, adjusted_FWHM, D,...
                             truncation, adjust_kernel );
        smooth_data = smooth_data / dx^2;
    elseif derivtype == 1 % Calculates the derivatives of the convolution field
        smooth_data(1,:,:,:) = fconv( expanded_lat_data,...
                                      { @(x) Gkerderiv( x, adjusted_FWHM ),...
                                        @(y) Gker( y, adjusted_FWHM ) },...
                                        D, truncation );
        smooth_data(2,:,:,:) = fconv( expanded_lat_data,...
                                      { @(x) Gker( x, adjusted_FWHM ),...
                                        @(y) Gkerderiv( y, adjusted_FWHM ) },...
                                        D, truncation );
        smooth_data = smooth_data / dx^3;
    else
        error( 'Higher derivatives are not supported' )  
    end
elseif D == 3
    % Increase the resolution of the raw data by introducing zeros
    expanded_lat_data = zeros( [ Dimhr, nsubj ] );
    expanded_lat_data( 1:( resAdd + 1 ):end, 1:( resAdd + 1 ):end, 1:( resAdd + 1 ):end, : ) = lat_data;
    
    % Adjust the FWHM to account for the dx
    adjusted_FWHM = FWHM/dx;
    
    % Run the smoothing
    if derivtype == 0 %Calculates the convolution field
        smooth_data = fconv( expanded_lat_data, adjusted_FWHM, D, truncation, adjust_kernel );
        smooth_data = smooth_data/dx^3;
    elseif derivtype == 1 %Calculates the derivatives of the convolution field
        smooth_data(1,:,:,:) = fconv(expanded_lat_data, {@(x)Gkerderiv(x,adjusted_FWHM), @(y)Gker(y,adjusted_FWHM), @(z)Gker(z,adjusted_FWHM)}, D, truncation);
        smooth_data(2,:,:,:) = fconv(expanded_lat_data, {@(x)Gker(x,adjusted_FWHM), @(y)Gkerderiv(y,adjusted_FWHM), @(z)Gker(z,adjusted_FWHM)}, D, truncation);
        smooth_data(3,:,:,:) = fconv(expanded_lat_data, {@(x)Gker(x,adjusted_FWHM), @(y)Gker(y,adjusted_FWHM), @(z)Gkerderiv(z,adjusted_FWHM)}, D, truncation);
        smooth_data = smooth_data/dx^3; % this needs to be checked! I think the scaling is wrong!
    else
        error('Higher derivatives are not supported')
    end
else
    error('D != 1,2,3 has not been implemented yet!')
end

end