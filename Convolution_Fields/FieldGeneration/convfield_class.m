function [ smooth_data, xvals_vecs, Kernel ] = convfield_class( lat_data,...
                                                         Kernel, resAdd, D,...
                                                         derivtype, enlarge )
% CONVFIELD_CLASS( lat_data, Kernel, resAdd, D, derivtype, enlarge )
% generates a convolution field from lattice data smoothed with an
% seperable kernel which can be specified. The generated field is evaluated
% on an equidistant grid with resolution increased by adding resAdd voxels
% inbetween each voxel in each dimension.
%
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  lat_data  data array T_1 x ... x T_D x N. Last index enumerates the
%            samples. Note that N > 1 is required!
%  Kernel    either an object of class SepKernel or a numeric.
%            If class SepKernel:
%              if derivtype = 0: the fields 'kernel' and 'truncation' must
%                                be specified.
%              if derivtype = 1: the fields 'dkernel' and 'dtruncation'
%                                must be specified.
%
%            If Kernel is numeric, the convolution field is generated by 
%            smoothing with an isotropic Gaussian kernel with FWHM = Kernel.
%            Truncation and adjust_kernel are set to be default values.
% Optional
%  resAdd     the amount of voxels added equidistantly inbetween the
%             existing voxels. Default is 1.
%  D          the dimension of the data, if this is left blank it is
%             assumed that nsubj = 1 and that the convolution field has 
%             the same numberof dimesions as lat_data
%  derivtype  0/1/2, 0 returns the convolution field, 1 its derivative
%             and 2 its second derivative (at all points). Default is 0
%             i.e to return the field. Note that if D > 1, then derivtype = 2
%             has not yet been implemented.
%  enlarge    a numeric which must be a positive integer or zero. The
%             convolution field is computed on a domain enlarged in each
%             direction by 'enlarge' voxels. Note if resAdd ~=0 the added
%             voxels are in high resolution. Default 0. 
%--------------------------------------------------------------------------
% %%% 1D
% %% Smoothing with increased resolution
% nvox = 10; xvals = 1:nvox; FWHM = 2;
% lat_data = normrnd(0,1,1,nvox); resadd = 10;
% lat_field = fconv(lat_data, FWHM);
% plot(xvals, lat_field, 'o-')
% hold on
% [convolution_field, xvals_fine] = convfield_class( lat_data, FWHM, resadd, 1);
% plot(xvals_fine{1},convolution_field)
% 
% %% Multiple subjects
% nsubj = 3; nvox = 100;
% lat_data = normrnd(0,1,nvox,nsubj);
% convolution_field = convfield_class( lat_data, FWHM, 0, 1 );
% plot(1:nvox,convolution_field)
% 
% %% 1D derivatives
% nvox = 100; resadd = 10; h = (1/(resadd+1)); 
% lat_data = normrnd(0,1,nvox,1);
% [convolution_field, xvals_fine] = convfield_class( lat_data, FWHM, resadd, 1);
% deriv1 = convfield_class( lat_data, FWHM, resadd, 1, 1 );
% deriv2 = diff(convolution_field)/h;
% plot(xvals_fine{1}, deriv1)
% hold on 
% plot(xvals_fine{1}(1:end-1), deriv2, '--')
% 
% % 1D derivative (multiple subjects) (Some minor differences here!)
% nvox = 100; FWHM = 3; xvals = 1:nvox; lat_data = normrnd(0,1,1,nvox);
% aderiv = @(x) applyconvfield( x, lat_data, @(y) GkerMVderiv(y, FWHM) );
% deriv_struct = convfield_class( lat_data, FWHM, resadd, 1, 1 );
% deriv_conv = convfield( lat_data, FWHM, resadd, 1, 1 );
% deriv_struct(1), deriv_conv(1), aderiv(1)
% 
% %% 2D
% Dim = [25,25];
% lat_data = normrnd(0,1,Dim);
% smooth_data = convfield_class( lat_data, FWHM, 0, 2);
% fine_data = convfield_class( lat_data, FWHM, 3, 2); %Convolution eval
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
% lat_data = normrnd(0,1,Dim); resadd = 1;
% derivfield = convfield_class( lat_data, FWHM, resadd, 2, 1);
% surf(reshape(derivfield(:,:,1), spacep(Dim,resadd)))
% title('2D 1st partial derivative of the convolution field')
% 
% %% Showing that the derivatives work
% Dim = [5,5]; lat_data = normrnd(0,1,Dim);
% point = [3,3]'; resadd = 100; h = 1/(1+resadd);
% 
% spaced_point = spacep(point,resadd);
% derivfield = convfield_class( lat_data, FWHM, resadd, 2, 1);
% 
% % convolution derivatives
% convfield_derivatives = squeeze(derivfield(spaced_point(1), spaced_point(2),:))
% 
% % Derivative using applyconvfield_class
% aderiv = @(x) applyconvfield( x, lat_data, @(y) GkerMVderiv(y, FWHM)  );
% acfield_derivatives = aderiv(point)
% 
% % Illustration on a fine lattice (not to be used in practice)
% smoothfield100 = convfield_class( lat_data, FWHM, resadd, 2);
% % smoothfield100 = convfield( lat_data, FWHM, resAdd, 2, 0);
% partialderiv_finelat(1) = (smoothfield100(spaced_point(1)+1, spaced_point(2)) - smoothfield100(spaced_point(1),spaced_point(2)))/h;
% partialderiv_finelat(2) = (smoothfield100(spaced_point(1), spaced_point(2) + 1) - smoothfield100(spaced_point(1),spaced_point(2)))/h;
% fine_lattice_derivatives = partialderiv_finelat'
% % note that it doesn't match perfectly because it's still a discrete
% % approximation, but that's why we want to use derivfield in the first
% % place!
% 
% % SPM (i.e. lattice) estimates of the derivative (quite off!)
% smoothfield_spm = convfield_class( lat_data, FWHM, 0, 2);
% spm_derivs(1) = (smoothfield_spm(point(1)+1, point(2)) - smoothfield_spm(point(1),point(2)));
% spm_derivs(2) = (smoothfield_spm(point(1), point(2) + 1) - smoothfield_spm(point(1),point(2)));
% spm_lattice_derivatives = spm_derivs'
% 
% %% 2D derivatives (multiple subjects)
% Dim = [5,5]; nsubj = 20;
% lat_data = normrnd(0,1,[Dim, nsubj]);
% derivfield = convfield_class( lat_data, FWHM, 1, 2, 1)
% 
% %% 3D
% %Compare to SPM
% Dim = [10,10,10]; FWHM = 3;
% lat_data = normrnd(0,1,Dim);
% subplot(1,2,1)
% spm_smooth_field = zeros(Dim); 
% spm_smooth(lat_data, spm_smooth_field, FWHM)
% surf(spm_smooth_field(:,:,5));
% title('Lattice Eval')
% subplot(1,2,2)
% cfield = convfield_class( lat_data, FWHM, 0, 3); %Convolution eval
% surf(cfield(:,:,5))
% title('Convolution Field Eval (no smoothing)')
% 
% %% Fine evaluation
% Dim = [10,10,10];
% resadd = 10; D = length(Dim); FWHM = 3; 
% slice = Dim(end)/2; spaced_slice = spacep(slice, resadd);
% lat_data = normrnd(0,1,Dim);
% cfield = convfield_class( lat_data, FWHM, resadd, D); %Convolution eval
% twoDcfieldslice = cfield(:,:,spaced_slice);
% zlimits = [min(twoDcfieldslice(:))-0.1, max(twoDcfieldslice(:))+0.1];
% 
% subplot(1,2,1)
% spm_smooth_field = zeros(Dim); 
% spm_smooth(lat_data, spm_smooth_field, FWHM)
% surf(spm_smooth_field(:,:,slice));
% title('Lattice Eval')
% zlim(zlimits)
% 
% subplot(1,2,2)
% surf(cfield(:,:,spaced_slice))
% title('Convolution Field Eval (Convn)')
% zlim(zlimits)
% 
% %% Compare to applyconvfield
% lat_data = normrnd(0,1,Dim);
% acfield = @(x) applyconvfield(x, lat_data, FWHM);
% Dim = [10,10,10];
% D = length(Dim); FWHM = 3; resadd = 0;
% cfield = convfield_class( lat_data, FWHM, resadd, D); 
% acfield([5,5,5]')
% cfield(5,5,5)
% acfield([1,1,10]')
% cfield(1,1,10)
% 
% %% % 3D derivatives (1 subject)
% Dim = [5,5,5]; D = length(Dim); FWHM = 3;
% lat_data = normrnd(0,1,Dim); resadd = 18;
% derivfield = convfield_class( lat_data, FWHM, 0, D, 1);
% aderiv = @(x) applyconvfield( x, lat_data, @(y) GkerMVderiv(y, FWHM)  );
% dfeval = squeeze(derivfield(3,3,3,:))
% aceval = aderiv([3,3,3]')
% spacing = 1/(1+resadd);
% spaced_point = spacep( [3,3,3]', resadd);
% cfield_fine = convfield_class( lat_data, FWHM, resadd, D);
% pointeval = cfield_fine(spaced_point(1),spaced_point(2),spaced_point(3));
% plusxeval = cfield_fine(spaced_point(1)+1,spaced_point(2),spaced_point(3));
% plusyeval = cfield_fine(spaced_point(1),spaced_point(2)+1,spaced_point(3));
% pluszeval = cfield_fine(spaced_point(1),spaced_point(2),spaced_point(3)+1);
% derivx = (plusxeval - pointeval)/spacing;
% derivy = (plusyeval - pointeval)/spacing;
% derivz = (pluszeval - pointeval)/spacing;
% 
% fine_lat_deriv = [derivx,derivy,derivz]'
% 
% %% 3D derivatives (Multiple subjects)
% Dim = [5,5,5]; D = length(Dim); FWHM = 3; nsubj = 2;
% lat_data = normrnd(0,1,[Dim,nsubj]); resadd = 0;
% derivfields = convfield_class( lat_data, FWHM, resadd, D, 1)
% aderiv = @(x) applyconvfield( x, lat_data(:,:,:,2), @(y) GkerMVderiv(y, FWHM)  )
% dfeval = squeeze(derivfields(3,3,3,2,:))
% aceval = aderiv([3,3,3]')
% 
% %% Adjusting the field (1D)
% nvox = 10; D = 1; FWHM = 2; lat_data = normrnd(0,1,1,nvox);
% kernel = SepKernel( D, FWHM ); kernel.adjust = 0.1; resadd = 10;
% [smoothfield, xvals_vecs] =  convfield_class( lat_data, FWHM, resadd, D);
% [adjust_field, xvals_vecs_adjust] = convfield_class( lat_data, kernel, 0, D);
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
% lat_data = normrnd(0,1,Dim); resadd = 9;
% [smoothfield, xvals_vecs] = convfield_class( lat_data, FWHM, resadd, D );
% 
% kernel = SepKernel( D, FWHM ); kernel.adjust = [0.1,0,0];
% [adjust_field, xvals_vecs_adjust] = convfield_class( lat_data, kernel, 0, D );
% 
% point = [1.1,1,1]'; spaced_point = spacep(point, resadd);
% plot(xvals_vecs{1}, smoothfield(:,spaced_point(2), spaced_point(3)))
% hold on
% plot(xvals_vecs_adjust{1},adjust_field(:,point(2),point(3)))
% acfield = @(tval) applyconvfield(tval, lat_data, FWHM);
% adjust_field(1,1)
% acfield(point)
% 
% %% %% Enlarging the field
% 
% %% 1D enlargement
% nvox = 10; lat_data = normrnd(0,1,1,nvox); FWHM = 3; D = 1; resadd = 0;
% cfield = convfield_class( lat_data', FWHM, resadd, D, 0, 1);
% acfield = @(tval) applyconvfield(tval, lat_data, FWHM);
% 
% cfield(1)
% acfield(0)
% 
% %% 1D enlargement
% nvox = 10; lat_data = normrnd(0,1,1,nvox); FWHM = 3; D = 1; resadd = 1;
% dx = 1/(1+resadd); enlarge = 1;
% cfield = convfield_class( lat_data', FWHM, resadd, D, 0, enlarge);
% acfield = @(tval) applyconvfield(tval, lat_data, FWHM);
% 
% cfield(1)
% acfield(1-dx*enlarge)
% %%
% nvox = 10; lat_data = normrnd(0,1,1,nvox); FWHM = 3; D = 1; resadd = 0;
% dx = 1/(1+resadd); enlarge = 5;
% cfield = convfield_class( lat_data', FWHM, resadd, D, 0, enlarge);
% acfield = @(tval) applyconvfield(tval, lat_data, FWHM);
% 
% cfield(1)
% acfield(1-dx*enlarge)
% 
% %% 2D enlargement
% Dim = [10,10]; lat_data = normrnd(0,1,Dim); FWHM = 3; D = 2; resadd = 1;
% dx = 1/(1+resadd); enlarge = 1;
% cfield = convfield_class( lat_data, FWHM, resadd, D, 0, enlarge)
% acfield = @(tval) applyconvfield(tval, lat_data, FWHM);
% 
% cfield(1,1)
% acfield([1-dx*enlarge,1-dx*enlarge]')
% 
% %% 3D enlargement
% Dim = [10,10,10]; lat_data = normrnd(0,1,Dim); FWHM = 3; D = 3; resadd = 1;
% dx = 1/(1+resadd); enlarge = 1;
% cfield = convfield_class( lat_data, FWHM, resadd, D, 0, enlarge);
% acfield = @(tval) applyconvfield(tval, lat_data, FWHM);
% 
% cfield(1,1,1)
% acfield([1-dx*enlarge,1-dx*enlarge,1-dx*enlarge]')
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport, Fabian Telschow                                              
%--------------------------------------------------------------------------


%% Get constants from the input lat_data field
%--------------------------------------------------------------------------

% Get size of the input data
slatdata = size( lat_data );
% Get number of dimensions of input data
D_latdata = length( slatdata );


%% Add/check optional values and define the kernel structure
%--------------------------------------------------------------------------

%%% D input
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

% Index for multidimensional coding
indexD = repmat( {':'}, 1, D );

%%% resAdd input
% Add resolution increase if missing
if ~exist( 'resAdd', 'var' )
    resAdd = 1;
end

% Get the difference between voxels with resolution increase
dx = 1 / ( resAdd + 1 );

% Reject input, if resAdd is to large in 3D
if D == 3 && ( resAdd > 18 )
    error( 'In 3D you shouldn''t use such high resolution for memory reasons' )
end

%%% derivtype input
if ~exist( 'derivtype', 'var' )
    derivtype = 0;
end

%%% enlarge input
if ~exist( 'enlarge', 'var' )
    enlarge = 0;
end

%% check and prepare the Kernel input structure
%--------------------------------------------------------------------------

if isnumeric( Kernel ) % If Kernel is numeric use an isotropic Gaussian Kernel
    
    % Set the FWHM parameter to be the user input adjusted by the resolution
    % factor
    FWHM = Kernel  / dx;
    
    % Change numerical Kernel input to structure
    Kernel = SepKernel( D, FWHM );

elseif isa( Kernel, 'SepKernel' )
    % Here probably a careful maniplualtion of the Kernel Class needs to
    % take place once we understand the factor.
else
    error( strcat( "The 'Kernel' must be either a numeric or an ",...
                   "object of class SepKernel!" ) );
end

%% Main function - separable case
%--------------------------------------------------------------------------

% Dimensions for domain of the field with increased resolution
Dimhr = ( Dim - 1 ) * resAdd + Dim; %Dimhr = Dim with high resolution

% Modify Dimhr by adding enlarge voxels to all sides
if enlarge ~= 0 
    if D == 1
        Dimhr = Dimhr + 2 * enlarge;
    else
        Dimhr = Dimhr + 2 * enlarge;
    end
end

% Setting up the default xvals_vecs
xvals_vecs  = cell( 1, D );

for d = 1:D
    xvals_vecs{d} = ( ( 1 - enlarge*dx ):dx:( Dim(d) + enlarge*dx ) )...
        + Kernel.adjust(d);
end

% Create index to fill the original data at the correct voxels of the high
% resolution data array
index = cell( [ 1 D ] );
for d = 1:D
    index{d} = ( enlarge + 1 ):( resAdd + 1 ):( Dimhr(d) - enlarge );
end

% Increase the resolution of the raw data by introducing zeros
expanded_lat_data = zeros( [ Dimhr, nsubj ] );
expanded_lat_data( index{:}, : ) = lat_data;

%%% Main loop: calculation of convolution fields
if D < 4
    % run the smoothing using fconv
    if derivtype == 0
        % calculates the convolution field
        smooth_data = fconv( expanded_lat_data, Kernel.kernel, D,...
                             Kernel.truncation(1), Kernel.adjust );
        % @Sam: can you explain the factor, naively I would have thought
        % that should be accounted for in the kernel and its derivative
        % (No I need to think about it but it makes it work lol, it's
        % something to do with taking FWHM = FWHM/dx)
        smooth_data = smooth_data / dx^D;
        
    elseif derivtype == 1
        % preallocate the output field for speed
        smooth_data = ones( [ Dimhr nsubj D ] );
        
        % get the gradient object of the Kernel
        dKernel = Gradient( Kernel );
        
        % calculates the derivatives of the convolution field
        for d = 1:D
            smooth_data(indexD{:},:,d) = fconv( expanded_lat_data,...
                                                dKernel.kernel{d}, D, ...
                                                dKernel.truncation(1) );
        end
        % @Sam: can you explain the factor 
        smooth_data = squeeze( smooth_data ) / dx^(D+1);
        
    else
        error( 'Higher derivatives are not supported' )  
    end
    
else
    error('D != 1,2,3 has not been implemented yet!')
end

%% Main function - non separable case
%--------------------------------------------------------------------------


return