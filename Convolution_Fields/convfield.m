function smooth_data = convfield( lat_data, FWHM, spacing, D, derivtype )
% CONVFIELD( lat_data, xvals_vecs, FWHM )
%--------------------------------------------------------------------------
% ARGUMENTS
% lat_data      the data on a lattice to be smoothed
% FWHM          the FWHM of the kernel with which to do smoothing
% spacing
% D             the dimension of the data
% derivtype     0/1/2, 0 returns the convolution field, 1 it's derivative
%               and 2 it's second derivative (at all points). Default is 0
%               i.e to return the field!
%--------------------------------------------------------------------------
% OUTPUT
%
%--------------------------------------------------------------------------
% EXAMPLES
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
% nsubj = 3;
% nvox = 100;
% lat_data = normrnd(0,1,nvox,nsubj);
% cfield = spm_conv(lat_data(:,1), FWHM);
% plot(1:nvox,cfield)
% hold on
% smooth_data = convfield( lat_data, FWHM, 1, 1 );
% plot(1:nvox,smooth_data(:,1))
% 
% %% 1D derivatives
% lat_data = normrnd(0,1,nvox,1);
% h = 0.01;
% smoothedfield = convfield( lat_data, FWHM, h, 1);
% deriv1 = convfield( lat_data, FWHM, h, 1, 1 );
% deriv2 = diff(smoothedfield)/h;
% plot(xvals_fine, deriv1 + 0.5)
% hold on 
% plot(xvals_fine(1:end-1), deriv2)
% 
% %% 2D
% Dim = [50,50];
% lat_data = normrnd(0,1,Dim)
% cfield = spm_conv(lat_data, FWHM)
% surf(cfield)
% smooth_data = convfield( lat_data, FWHM, 0, 0, 1)
% surf(smooth_data)
% 
% % 2D derivatives
% derivfield = convfield( lat_data, FWHM, 0, 1, 1)
% surf(reshape(derivfield(1,:), Dim))
% resAdd = 100;
% smoothfield100 = convfield( lat_data, FWHM, resAdd, 0, 1);
% derivfield100 = convfield( lat_data, FWHM, resAdd, 1, 1);
% point = [500,500]
% ((smoothfield100(point(1), point(2) + 1) - smoothfield100(point(1),point(2)))/(1/(resAdd +1)))
% ((smoothfield100(point(1)+1, point(2)) - smoothfield100(point(1),point(2)))/(1/(resAdd +1)))
% derivfield100(:,point(1), point(2))
% % note that it's still not perfect because it's still a discrete
% % approximation, but that's why we want to use derivfield in the first
% % place!!
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
if nargin < 3
    spacing = 0.1;
end
if nargin < 5
    derivtype = 0;
end

slatdata = size(lat_data);
D_latdata = length(slatdata);
if D > 1
    if D_latdata == D
        nsubj = 1;
        Dim = slatdata;
    else
        nsubj = slatdata(end);
        Dim = slatdata( 1 : end-1 );
    end
else
    vert2horz = 0;
    if D_latdata == 2 && slatdata(1) == 1
        nsubj = 1;
        Dim = slatdata(slatdata > 1);
        vert2horz = 1; % I.e. if a horizontal field is entered it returns one as well
    elseif D_latdata == 2 && slatdata(2) == 1
        nsubj = 1;
        Dim = slatdata(slatdata > 1);
    else %I.e. if D_lat_data == 2 and slatdata(1) and slatdata(2) are both > 1
        % or if D_lat_data > 2 and there are some 1 dimensional dimensions
        % hanging around for some reason.
        slatdata_squeezed = size(squeeze(lat_data));
        Dim = slatdata_squeezed(1);
        nsubj = slatdata_squeezed(2);
    end
end

resAdd = floor(1/spacing-1);
dx = 1/(resAdd+1); %Gives the difference between voxels

% Dimensions for field with increased resolution
Dimhr = ( Dim - 1 ) * resAdd + Dim; %Dimhr = Dim with high resolution.

truncation = ceil( 4*FWHM2sigma(FWHM) );
gridside  = -truncation:dx:truncation;

% convolution kernel and derivatives to be used with convn
if D == 1
    % increase the resolution of the raw data by introducing zeros
    expanded_lat_data = zeros( [ Dimhr, nsubj] );
    expanded_lat_data( 1:(resAdd + 1):end, : ) = lat_data;
    
    % grid for convolution kernel
    
    if derivtype == 0
        h = Gker( gridside, FWHM, 1 );
        smooth_data = convn( expanded_lat_data', h, 'same' )';
    elseif derivtype == 1
        [ ~, dxh ] = Gker( gridside, FWHM, 1 );
        smooth_data = convn( expanded_lat_data', dxh, 'same' )';
    elseif derivtype == 2
        [ ~, ~, d2xh ] = Gker( gridside, FWHM, 1 );
        smooth_data = convn( expanded_lat_data', d2xh, 'same' )';
    else
        error('Higher derivatives are not supported')
    end
    
    if vert2horz && nsubj == 1
        smooth_data = smooth_data';
    end
elseif D == 2
    expanded_lat_data = zeros( [ Dimhr, nsubj ] );
    expanded_lat_data( 1:( resAdd + 1 ):end, 1:( resAdd + 1 ):end, : ) = lat_data;
    
    % grid for convolution kernel
    [x,y] = meshgrid( gridside, gridside );
    xvals = [x(:), y(:)]';
    
    % convolution kernels to be used with convn
    if derivtype == 0
        h   = reshape( GkerMV( xvals, FWHM ), size(x) );
        smooth_data  = convn( expanded_lat_data, h, 'same' );
    elseif derivtype == 1
        smooth_data = zeros( [D, size(expanded_lat_data)]);
        dh  = GkerMVderiv( xvals, FWHM );
        dxh = reshape( dh(1,:), size(x) );
        dyh = reshape( dh(2,:), size(x) );
        smooth_data(1,:,:)  = convn( expanded_lat_data, dxh, 'same' );
        smooth_data(2,:,:) = convn( expanded_lat_data, dyh, 'same' );
    else
        error('Higher derivatives are not supported')  
    end
end

end

