function smooth_data = cfield( lat_data, FWHM, resAdd, derivtype, nsubj )
% CFIELD( lat_data, xvals_vecs, FWHM )
%--------------------------------------------------------------------------
% ARGUMENTS
% lat_data
% FWHM
% resAdd     the number of points between each voxel. Default is zero
%            i.e. to just compute the smooth field on the lattice.
%--------------------------------------------------------------------------
% OUTPUT
%
%--------------------------------------------------------------------------
% EXAMPLES
% %%% 1D
% % Smoothing with increased resolution
% nvox = 100; xvals = 1:nvox;
% FWHM = 3;
% lat_data = normrnd(0,1,1,nvox)
% convfield = inter_conv1D( lat_data, FWHM, 0.01);
% plot(1:0.01:nvox,convfield)
% hold on
% smooth_data = cfield( lat_data', 99, FWHM )
% plot(1:0.01:nvox,smooth_data)
%
% %Smoothing with the same resolution
% lat_data = normrnd(0,1,1,nvox)
% convfield = spm_conv(lat_data, FWHM)
% plot(1:nvox,convfield)
% hold on
% smooth_data = cfield( lat_data', FWHM )
% plot(1:nvox, smooth_data)
%
% % Multiple subjects
% nsubj = 3;
% lat_data = normrnd(0,1,nvox,nsubj)
% convfield = spm_conv(lat_data(:,1), FWHM)
% plot(1:nvox,convfield)
% hold on
% smooth_data = cfield( lat_data, FWHM )
% plot(1:nvox,smooth_data(:,1))
%
% 1D derivatives
% deriv1 = cfield( lat_data(:,1), FWHM, resAdd, 1, nsubj );
% deriv2 = diff(smoothedfield(:,1));
% plot(1:nvox, deriv1)
% hold on 
% plot(1:(nvox-1), deriv2)
%
% %%% 2D
% Dim = [50,50];
% lat_data = normrnd(0,1,Dim)
% convfield = spm_conv(lat_data, FWHM)
% surf(convfield)
% smooth_data = cfield( lat_data, FWHM, 0, 0, 1)
% surf(smooth_data)
%
% % 2D derivatives
% derivfield = cfield( lat_data, FWHM, 0, 1, 1)
% surf(reshape(derivfield(1,:), Dim))
% resAdd = 100;
% smoothfield100 = cfield( lat_data, FWHM, resAdd, 0, 1);
% derivfield100 = cfield( lat_data, FWHM, resAdd, 1, 1);
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
    resAdd = 0;
end
if nargin < 4
    derivtype = 0;
end

if nargin >= 5 && nsubj == 1
    Dim = size(lat_data);
else
    sY = size(lat_data);
    nsubj = sY(end);
    Dim = sY( 1 : end-1 );
end
D = length(Dim);

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
elseif D == 2
    expanded_lat_data = zeros( [ Dimhr, nsubj ] );
    expanded_lat_data( 1:( resAdd + 1 ):end, 1:( resAdd + 1 ):end, : ) = lat_data;
    
    % grid for convolution kernel
    [x,y] = meshgrid( gridside, gridside );
    xvals = [x(:), y(:)]';
    
    % convolution kernels to be used ith convn
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

