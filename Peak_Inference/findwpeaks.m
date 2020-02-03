function peak_locs = findwpeaks(lat_data, Kernel, xvals_vecs, peak_est_locs, mask)
% FINDWPEAKS( lat_data, Kprime, xvals_vecs, peak_est_locs, Kprime2, truncation, mask )
% calculates the locations of peaks in a convolution field using Newton Raphson.
%--------------------------------------------------------------------------
% ARGUMENTS
% lat_data      A D by nsubj matrix array giving the value of the lattice 
%               field at every point.
% Kprime        a function handle giving the derivative of the kernel. If
%               this is numeric the kernel is taken to be Gaussian with the 
%               value as the FWHM. I.e. if Kprime = 2, then within the 
%               function Kprime is set to be @(x) GkerMVderiv(x,2).
% xvals_vecs    a D-dimensional cell array whose entries are vectors giving the
%               xvalues at each each dimension along the lattice. It assumes
%               a regular, rectangular lattice (though within a given
%               dimension the voxels can be spaced irregularly).
%               I.e suppose that your initial lattice grid is a
%               4by5 2D grid with 4 voxels in the x direction and 5 in
%               the y direction. And that the x-values take the values:
%               [1,2,3,4] and the y-values take the values: [0,2,4,6,8].
%               Then you would take xvals_vecs = {[1,2,3,4], [0,2,4,6,8]}.
%               The default is to assume that the spacing between the
%               voxels is 1. If only one xval_vec direction is set the
%               others are taken to range up from 1 with increment given by
%               the set direction.
% peak_est_locs a D by npeaks matrix giving the initial estimates of the 
%               location of the peaks. If this is instead an integer: top  
%               then the top number of maxima are considered and initial 
%               locations are estimated from the underlying data. If this 
%               is not specified then it is set to 1, i.e. only considering
%               the maximum.
% Kprime2       the 2nd derivative of the kernel. If this is not set then
%               if the initial Kernel is Gaussian it is taken to be the 2nd
%               derivative of the corresponding Gaussian kernel. Otherwise
%               it is estimated from Kprime.
% truncation
% mask
%--------------------------------------------------------------------------
% OUTPUT
% peak_locs   the true locations of the top peaks in the convolution field.
%--------------------------------------------------------------------------
% EXAMPLES
% % 1D
% L = 100;
% nsubj = 50;
% lat_data = normrnd(0,1, nsubj, L)';
% findtpeaks(lat_data, 3)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport.
Ldim = size(lat_data);
Ldim = Ldim(1:end-1);
D = length(Ldim);

if nargin < 5
    if D == 1
        mask = ones(1,Ldim);
    else
        mask = ones(Ldim);
    end
end

if isnan(sum(lat_data(:)))
   error('Cant yet deal with nans') 
end

%Setting up xvals_vecs
if nargin < 4
    xvals_vecs = {1:Ldim(1)};  %The other dimensions are taken case of below.
end
if ~iscell(xvals_vecs)
    xvals_vecs = {xvals_vecs};
end
if length(xvals_vecs) < D
    increm = xvals_vecs{1}(2) - xvals_vecs{1}(1);
    for d = (length(xvals_vecs)+1):D
        xvals_vecs{d} = xvals_vecs{1}(1):increm:(xvals_vecs{1}(1)+increm*(Ldim(d)-1));
    end
end

% nsubj = size(lat_data, 1);
xvals_vecs_dims = zeros(1,D);
for d = 1:D
    xvals_vecs_dims(d) = length(xvals_vecs{d});
end
if ~isequal(xvals_vecs_dims, Ldim)
    error('The dimensions of xvals_vecs must match the dimensions of lat_data')
end

wcf = @(tval) applyconvfield( tval, lat_data, xvals_vecs, Kernel )./applyconvfield( tval, mask, xvals_vecs, Kernel );
h = 0.0001;
if D == 1
    wcf_deriv = @(tval) (wcf(tval+h) - wcf(tval))/h;
    wcf_deriv2 = @(tval) (wcf_deriv(tval+h) - wcf_deriv(tval))/h;
elseif D == 2
    wcf_deriv = @(x) [(wcf(x+h*[1,0]') - wcf(x))/h, (wcf(x+h*[0,1]') - wcf(x))/h]';
    wcf_deriv2 = @(x) [(wcf_deriv(x+h*[1,0]') - wcf_deriv(x))/h, (wcf_deriv(x+h*[0,1]') - wcf_deriv(x))/h];
elseif D == 3
    wcf_deriv = @(x) [(wcf(x+h*[1,0,0]') - wcf(x))/h, (wcf(x+h*[0,1,0]') - wcf(x))/h, (wcf(x+h*[0,0,1]') - wcf(x))/h]';
    wcf_deriv2 = @(x) [(wcf_deriv(x+h*[1,0,0]') - wcf_deriv(x))/h, (wcf_deriv(x+h*[0,1,0]') - wcf_deriv(x))/h, (wcf_deriv(x+h*[0,0,1]') - wcf_deriv(x))/h];
end

if nargin < 4 
    peak_est_locs = 1; %I.e. just consider the maximum.
end

% At the moment this is just done on the initial lattice. Really need to
% change so that it's on the field evaluated on the lattice.

if isequal(size(peak_est_locs), [1,1])
    top = peak_est_locs;
    xvalues_at_voxels = xvals2voxels( xvals_vecs );
    if D < 3
        smoothed_data = wcf(xvalues_at_voxels);
        smoothed_data = reshape(smoothed_data, size(mask));
        %Using the above no longer need the Ktype condition
        %             smoothed_data = spm_conv(lat_data, FWHM);
    elseif D == 3
        smoothed_data = zeros(size(lat_data));
        weights = zeros(size(lat_data));
        spm_smooth(lat_data, smoothed_data, FWHM)
        spm_smooth(mask, weights, FWHM);
        smoothed_data = smoothed_data./weights;
    else
        error('Not yet ready for 4D or larger images')
    end
    max_indices = lmindices(smoothed_data, top, mask)'; %Note the transpose here! It's necessary for the input to other functions.
    if D == 1
        max_indices = max_indices';
    end
    top = length(max_indices);  
    peak_est_locs = zeros(D, top);
    for I = 1:D
        peak_est_locs(I, :) = xvals_vecs{I}(max_indices(I,:));
    end
end
if D == 1 && isnan(peak_est_locs(1))
    peak_est_locs = peak_est_locs(2:end);
end
npeaks = size(peak_est_locs, 2);

peak_locs = zeros(D, npeaks);
for peakI = 1:npeaks
    tol = min(abs(wcf_deriv(peak_est_locs(:,peakI)))/100000, 0.0001);
    notconverged = 1;
    while notconverged > 0 %This loop is to deal with situations where things don't converge on the first initialization.
        try
            peak_locs(:, peakI) = NewtonRaphson(wcf_deriv, peak_est_locs(:, peakI) + (notconverged-1)*0.1, wcf_deriv2, tol);
%             peak_locs(peakI) = NewtonRaphson(fderiv, peak_est_locs(peakI) + (notconverged-1)*0.1, fderiv2, tol);
            notconverged = 0;
        catch
            notconverged = notconverged + 1;
            if notconverged > 10
                error('notconverged has not converged')
            end
        end
    end
end

end