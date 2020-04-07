function peak_locs = findconvpeaks(lat_data, Kernel, peak_est_locs, mask, xvals_vecs, truncation)
% FINDCONVPEAKS( lat_data, Kprime, xvals_vecs, peak_est_locs, Kprime2, truncation, mask )
% calculates the locations of peaks in a convolution field using Newton Raphson.
%--------------------------------------------------------------------------
% ARGUMENTS
% lat_data      A D-dimensional matrix array giving the value of the lattice 
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
%               the maximum. If D = 1 and you wish to specify multiple
%               peaks rather than a number of peaks then you need to begin
%               your input as [NaN, peakestloc1, peakestloc2, ...].
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
% These examples start on negative definite areas so it is easy to then use
% Newton-Raphson to find a maximum. It's more difficult if that's not the
% case!!
% % 1D 
% Y = [1,2,1];
% findconvpeaks(Y, 3, 1)
%
% % 2D
% Y = [1,1,1,1;1,2,2,1;1,2,2,1;1,1,1,1];
% findconvpeaks(Y, 3, [2,2]')
%
% Y = [1,1,1;2,1,1;1,1,1] %I.e so the peak will be outside the mask!
% mask = [0,1,1;0,1,1;0,1,1];
% findconvpeaks(Y, 3, [2,2]', mask)
%
% Y = [1,1,1;10,1,1;1,1,1] %I.e so the peak will be outside the mask!
% mask = [0,1,1;0,1,1;0,1,1];
% findconvpeaks(Y, 3, [3,3]', mask) %Need to fix so that
% initialization at [2,2] is fine as well!
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport.
Ldim = size(lat_data);
D = length(Ldim);
if Ldim(2) == 1
    lat_data = lat_data';
    Ldim = Ldim';
end
if Ldim(1) == 1
    D = D - 1;
    Ldim = Ldim(2:end);
    if D > 1
        lat_data = reshape(lat_data, Ldim);
    end
end
if size(peak_est_locs, 1) ~= D
    error('peak_est_locs is the wrong dimension!')
end
if nargin < 3 
    peak_est_locs = 1; %I.e. just consider the maximum.
end
if nargin < 4
    if D == 1
        mask = ones(1,Ldim);
    else
        mask = ones(Ldim);
    end
end
if nargin < 6
    truncation = 0;
end

if isnan(sum(lat_data(:)))
   error('Cant yet deal with nans') 
end

Ktype = '';  
if isnumeric(Kernel)
    FWHM = Kernel;
    Ktype = 'G';
end

%Setting up xvals_vecs
if nargin < 5
    xvals_vecs = {1:Ldim(1)}; %The other dimensions are taken case of below.
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

field = @(tval) applyconvfield(tval, lat_data, FWHM, truncation, xvals_vecs );

% At the moment this is just done on the initial lattice. Really need to
% change so that it's on the field evaluated on the lattice.

if isequal(size(peak_est_locs), [1,1]) && floor(peak_est_locs(1)) == peak_est_locs(1)
    top = peak_est_locs;
    if strcmp(Ktype, 'G')
        if D < 3
            xvalues_at_voxels = xvals2voxels(xvals_vecs);
            smoothed_data = applyconvfield(xvalues_at_voxels, lat_data, FWHM, truncation, xvals_vecs);
            smoothed_data = reshape(smoothed_data, size(mask));
            %Using the above no longer need the Ktype condition
%             smoothed_data = spm_conv(lat_data, FWHM);
        elseif D == 3
            smoothed_data = zeros(size(lat_data));
            spm_smooth(lat_data, smoothed_data, FWHM);
        else
            error('Not yet ready for 4D or larger images')
        end
        [~,~,max_indices] = lmindices(smoothed_data, top, mask); %Note the transpose here! It's necessary for the input to other functions.
        max_indices = max_indices';
    else
        [~,~,max_indices] = lmindices(lat_data, top, mask); %In 3D may need to use spm_smooth for this bit!
        max_indices = max_indices';
    end
    if D == 1
        max_indices = max_indices';
    end
    top = size(max_indices, 2);  
    peak_est_locs = zeros(D, top);
    for d = 1:D
        peak_est_locs(d, :) = xvals_vecs{d}(max_indices(d,:));
    end
end
if D == 1 && isnan(peak_est_locs(1)) %This allows for multiple 1D peaks!
    peak_est_locs = peak_est_locs(2:end);
end
npeaks = size(peak_est_locs, 2);

peak_locs = zeros(D, npeaks);
A = [eye(D);-eye(D)];
% b = [Ldim(:)+0.5;ones(D,1)-0.5];
b = zeros(2*D,1);
for d = 1:D
    b(d) = xvals_vecs{d}(end);
end
for d = 1:D
    b(d+D) = -xvals_vecs{d}(1);
end
%Need to discuss which boundary to use with Fabian!!!

options = optimoptions(@fmincon,'Display','off'); %Ensures that no output is displayed.
for peakI = 1:npeaks
    peak_locs(:, peakI) = fmincon(@(tval) -field(tval), peak_est_locs(:, peakI), A, b, [], [], [], [], [], options);
end

end


% for peakI = 1:npeaks
%     %     applyconvfield_gen(peak_est_locs(:, peakI), lat_data, Kprime, xvals_vecs )
%     %     field_deriv(peak_est_locs(:, peakI))
%     %     tol = min(abs(field_deriv(peak_est_locs(:,peakI)))/100000, 0.0001);
%     if nargin < 6
%         tol = min(abs(field_deriv(peak_est_locs(:,peakI)))/100000, 0.0001);
%     end
%     peak_locs(:, peakI) = findpeak(peak_est_locs(:, peakI), field_deriv, field_deriv2, mask, 1, tol, 0.05, 0.01, field);
%     if field(peak_locs(:, peakI)) < field(peak_est_locs(:, peakI))
%         peak_locs(:, peakI) = findpeak(peak_est_locs(:, peakI), field_deriv, field_deriv2, mask, 1, tol, 0.01, 0.0001, field);
%     end
%     % %     peak_locs(:, peakI) = gascent( peak_est_locs(:, peakI), field_deriv, 0.01, tol, field);
% %     peak_locs(:, peakI) = NewtonRaphson(field_deriv, peak_est_locs(:, peakI), field_deriv2, tol);
% end

% peak_locs(:, peakI) = findpeak( peak_est_locs(:, peakI), field_deriv, fprime, fprime2, 1, tol );
% %     peak_locs(:, peakI) = NewtonRaphson(field_deriv, peak_est_locs(:, peakI), field_deriv2, tol);
%     if (isnan(sum(peak_locs(:, peakI))) || norm(peak_locs(:, peakI) - peak_est_locs(:, peakI)) > 3)&& strcmp(Ktype, 'G') 
%         ninter = 0.25; %Could be adjusted starting bigger and made to get smaller?
%         subset_xvals_vecs = cell(1,D);
%         for d = 1:D
%             subset_xvals_vecs{d} = (peak_est_locs(d, peakI) - 1 + ninter):ninter:(peak_est_locs(d, peakI) + 1 - ninter);
%         end
%         subset_xvaluesatvoxels = xvals2voxels(subset_xvals_vecs);
%         field_around_peak = field(subset_xvaluesatvoxels);
%         [~,max_index] = max(field_around_peak);
%         new_est_peak_loc = subset_xvaluesatvoxels(:, max_index);
%         peak_locs(:, peakI) = NewtonRaphson(field_deriv, new_est_peak_loc, field_deriv2, tol);
%     end