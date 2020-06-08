function [peaklocs, peakvals] = findconvpeaks(lat_data, Kernel, peak_est_locs, mask, xvals_vecs, truncation)
% FINDCONVPEAKS(lat_data, Kernel, peak_est_locs, mask, xvals_vecs, truncation)
% calculates the locations of peaks in a convolution field.
%--------------------------------------------------------------------------
% ARGUMENTS
% lat_data      A D-dimensional matrix array giving the value of the lattice 
%               field at every point.
% Kernel        the smoothing Kernel (as a function). If Kernel is a postive
%               number then an isotropic Gaussian kernel with dimension D and
%               FWHM = Kernel is used
% xvals_vecs    a D-dimensional cell array whose entries are vectors giving the
%               xvalues at each each dimension along the lattice. It assumes
%               a regular, rectangular lattice (though within a given
%               dimension the voxels can be spaced irregularly).
%               I.e suppose that your initial lattice grid is a
%               4by5 2D grid with 4 voxels in the x direction and 5 in
%               the y direction the y-values take the values: [0,2,4,6,8].
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
% % 1D 
% Y = [1,2,1];
% findconvpeaks(Y, 3, 1)
% 
% % 1D with different xvals_vecs
% Y = [1,1,2,2,1,1];
% xvals_vecs = 11:(length(Y)+10);
% findconvpeaks(Y, 3, 1, ones(1,length(Y)), xvals_vecs)
% 
% % 1D multiple peaks
% Y = [1,2,1,1,1,1,2,1];
% findconvpeaks(Y, 3, 2) % Top 2 peaks
% findconvpeaks(Y, 3, [NaN, 1,6]) % Top 2 peaks starting initializing at 1, 6, 
%                                   needs the NaN in 1D to differentiate from 
%                                   the top number of peaks 
% findconvpeaks(Y, 3, [NaN, 4.5]) returns 4.5 which is not the max so watch
%                                   out it can get stuck atm!
% 
% % 2D
% Y = [1,1,1,1;1,2,2,1;1,2,2,1;1,1,1,1];
% [maxloc, maxval] = findconvpeaks(Y, 2, 1)
% fine_eval = convfield(Y, 2, 0.01, 2);
% max(fine_eval(:))
% 
% % With masking!
% FWHM = 3;
% Y = [10,1,1;1,1,1;10,1,1]; %I.e. so the peak will be outside the mask!
% mask = [1,1,1;0,1,1;1,1,1];
% findconvpeaks(Y, FWHM, [2,2]', mask)
% field = @(tval) applyconvfield(tval, Y, FWHM, -1, xvals_vecs, mask);
% 
% % View masked field:
% masked_field = @(x) mask_field(x, mask).*field(x);
% fieldeval = reshape(masked_field(xvaluesatvoxels), [length(xvals_fine),length(xvals_fine)]);
% surf(fieldeval)
% 
% Y = [1,1,1;10,1,1;1,1,1] %I.e so the peak will be outside the mask!
% mask = [0,1,1;0,1,1;0,1,1];
% findconvpeaks(Y, 3, [3,3]', mask) %Need to fix so that
% initialization at [2,2] is fine as well!
% 
% %2D multiple peaks
% Y = [5,1,1,1;1,1,1,1;1,1,1,1;1,1,1,5]
% surf(convfield(Y, 2, 0.1, 2))
% findconvpeaks(Y, 2, [1,1;4,4]')
%
% %% The maximum can lie on the boundary (could set this as an exercise)
% FWHM = 3; xvals_fine = 1:0.1:3;
% xvaluesatvoxels = xvals2voxels({xvals_fine,xvals_fine});
% Y = [10,1,1;1,1,1;10,1,1]; %I.e. so the peak will be outside the mask!
% mask = [1,1,1;0,1,1;1,1,1];
% field = @(tval) applyconvfield(tval, Y, FWHM, -1, xvals_vecs, mask);
% mfield = @(x) mask_field(x, mask);
% masked_field = @(x) mfield(x).*field(x);
% fieldeval = reshape(masked_field(xvaluesatvoxels), [length(xvals_fine),length(xvals_fine)]);
% surf(fieldeval)
% findconvpeaks(Y, FWHM, [2,2]', mask)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
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
if size(peak_est_locs, 1) ~= D && ~isequal(size(peak_est_locs), [1,1])
    error('peak_est_locs is the wrong dimension!')
end
if nargin < 3 
    peak_est_locs = 1; %I.e. just consider the maximum.
end
if ~exist('mask', 'var') || isequal(mask, NaN)
    if D == 1
        mask = ones(1,Ldim);
    else
        mask = ones(Ldim);
    end
end
if nargin < 6
    truncation = -1;
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
if nargin < 5 || isequal(xvals_vecs, NaN)
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

xvals_vecs_dims = zeros(1,D);
for d = 1:D
    xvals_vecs_dims(d) = length(xvals_vecs{d});
end
if ~isequal(xvals_vecs_dims, Ldim)
    error('The dimensions of xvals_vecs must match the dimensions of lat_data')
end

field = @(tval) applyconvfield(tval, lat_data, FWHM, truncation, xvals_vecs, mask );
if ~isequal(mask, ones(Ldim)) && ~isequal(mask, ones([1, Ldim]))
    masked_field = @(x) mask_field( x, mask, xvals_vecs ).*field(x);
else
    masked_field = field;
end

% warning('need masking here! and for the maximization!! though mahel')

% At the moment this is just done on the initial lattice. Really need to
% change so that it's on the field evaluated on the lattice.

if isequal(size(peak_est_locs), [1,1]) && floor(peak_est_locs(1)) == peak_est_locs(1)
    top = peak_est_locs;
    if strcmp(Ktype, 'G')
        if D < 4
            smoothed_data = fconv( lat_data, FWHM , D );
        else
            error('Not yet ready for 4D or larger images')
        end
        max_indices = lmindices(smoothed_data, top, mask);
    else
        max_indices = lmindices(lat_data, top, mask); %In 3D may need to use spm_smooth for this bit!
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

% Obtain the boundary of the mask
boundary = bndry_voxels( logical(mask), 'full' );

% Obtain the box sizes within which to search for the maximum
% Assign boxsize of 0.5 for voxels on the boundary and 1.5 for voxels not
% on the boundary.
npeaks = size(peak_est_locs, 2); % Calculate the number of estimates
s_mask = size(mask);             % Obtain the size of the mask
box_sizes = 1.5*ones(1,npeaks);
for I = 1:npeaks
    converted_index = convind( peak_est_locs(:,I), s_mask );
    if boundary(converted_index)
        box_sizes(I) = 0.5;
    end
end

% Find local maxima
[ peaklocs, peakvals ] = findlms( masked_field, peak_est_locs, box_sizes );

end


% DEPRECATED
% npeaks = size(peak_est_locs, 2);
% peakvals = zeros(1,npeaks);
% peaklocs = zeros(D, npeaks);
% 
% A = [eye(D);-eye(D)];
% % b = [Ldim(:)+0.5;ones(D,1)-0.5];
% b = zeros(2*D,1);
% for d = 1:D
%     b(d) = xvals_vecs{d}(end);
% end
% for d = 1:D
%     b(d+D) = -xvals_vecs{d}(1);
% end
% %Need to discuss which boundary to use with Fabian!!!
% 
% options = optimoptions(@fmincon,'Display','off'); %Ensures that no output is displayed.
% for peakI = 1:npeaks
%     peaklocs(:, peakI) = fmincon(@(tval) -field(tval), peak_est_locs(:, peakI), A, b, [], [], [], [], [], options);
%     peakvals(peakI) = field(peaklocs(:, peakI));
% end

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