function [peaklocs, peakvals] = findconvpeaks_t(lat_data, Kernel, peak_est_locs, mask, xvals_vecs, truncation )
% FINDCONVPEAKS_T( lat_data, Kprime, xvals_vecs, peak_est_locs, Kprime2, truncation, mask )
% calculates the locations of peaks in a convolution field using Newton Raphson.
%--------------------------------------------------------------------------
% ARGUMENTS
% lat_data      A D by nsubj matrix array giving the value of the lattice
%               field at every point.
% Kernel
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
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Calculate the size of the data including number of subjects
Ldim = size(lat_data);

% Calculate the dimensions of the data (as Ldim(end) = nsubj)
Ldim = Ldim(1:end-1); % Since for t-fields we always assume at least 2 subjects
D = length(Ldim); 

%%  add/check optional values
%--------------------------------------------------------------------------
if ~exist('peak_est_locs', 'var')
    peak_est_locs = 1; %Default is to just consider the maximum.
end
if ~exist('mask', 'var') || isequal(mask, NaN)
    if D == 1
        mask = ones(Ldim,1);
%         mask = ones(1,Ldim);
    else
        mask = ones(Ldim);
    end
end
if ~exist('truncation', 'var')
    truncation = -1;
end

% Setting up xvals_vecs
if ~exist('xvals_vecs', 'var')
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

xvals_vecs_dims = zeros(1,D);
for d = 1:D
    xvals_vecs_dims(d) = length(xvals_vecs{d});
end

%% Error checking
%--------------------------------------------------------------------------
if isnan(sum(lat_data(:)))
    error('Can''t yet deal with nans')
end
if ~isequal(xvals_vecs_dims, Ldim)
    error('The dimensions of xvals_vecs must match the dimensions of lat_data')
end

%%  main function
%--------------------------------------------------------------------------
% Define t convolution field
tcf = @(tval) applyconvfield_t( tval, lat_data, Kernel, truncation, ...
                                                        xvals_vecs, mask );
% Mark the field (if necessary)
if ~isequal(mask, ones(Ldim)) && ~isequal(mask, ones([1, Ldim]))
    masked_field = @(x) mask_field( x, mask, xvals_vecs ).*tcf(x);
else
    masked_field = tcf;
end

% At the moment this is just done on the initial lattice. Need to
% change so that it's on the field evaluated on the lattice.
if (D > 2 || isnumeric(peak_est_locs))  && (isequal(size(peak_est_locs),...
                   [1,1])) && (floor(peak_est_locs(1)) == peak_est_locs(1))
    % Set the number of maxima to define
    numberofmaxima2find = peak_est_locs;
    
    % Calculate the t field on the lattice
    teval_lat = convfield_t( lat_data, Kernel, 0 ); 
    
    % Find the top local maxima of the field on the lattice
    max_indices = lmindices(teval_lat, numberofmaxima2find, mask); 
    
    % In D = 1 you need to transpose (CHECK THIS) Note the transpose here! 
    % It's necessary for the input to other functions.
    if D == 1
        max_indices = max_indices'; 
    end
    
    % Reset this in case there are less maxima than desired
    numberofmaxima2find = size(max_indices,2);
    
    % Initialize a matrix to store peak locations
    peak_est_locs = zeros(D, numberofmaxima2find);
    for I = 1:D
        peak_est_locs(I, :) = xvals_vecs{I}(max_indices(I,:));
    end
elseif iscell(peak_est_locs)
    % For D = 1, if you want to calculate multiple peaks you need to enter
    % a cell array (to differentiate between this and the top number of
    % peaks).
    peak_est_locs = cell2mat(peak_est_locs);
end

% Obtain the boundary of the mask
boundary = bndry_voxels( logical(mask), 'full' );

% Obtain the box sizes within which to search for the maximum
% Assign boxsize of 0.5 for voxels on the boundary and 1.5 for voxels not
% on the boundary.
npeaks = size(peak_est_locs, 2); % Calculate the number of estimates
s_mask = size(mask);             % Obtain the size of the mask

% Set the default box sizes within which to search for maxima
box_sizes = repmat({1.5}, 1, npeaks);

% For peaks initialized on the boundary change their box sizes to 0.5
for I = 1:npeaks
    if D > 1
        converted_index = convind( peak_est_locs(:,I) - xvals_starts_at' + 1, s_mask );
    else
        converted_index = peak_est_locs(:,I) - xvals_starts_at' + 1;
    end
    if boundary(round(converted_index)) % Need to come back to this an make it more general
        box_sizes{I} = 0.5;
    end
end

% Find local maxima
[ peaklocs, peakvals ] = findlms( masked_field, peak_est_locs, box_sizes );

end

% DEPRECATED:
%     top = peak_est_locs;
% %     xvalues_at_voxels = xvals2voxels( xvals_vecs );
%     teval_lat = smoothtstat(lat_data, Kernel);
% %     if D < 3
% %         teval_lat = tcf(xvalues_at_voxels);
% %     else
% %         smoothed_field = zeros([Ldim, nsubj]);
% %         smoothing_store = zeros(Ldim);
% %         for subj = 1:nsubj
% %             spm_smooth(lat_data(:,:,:,subj), smoothing_store, Kernel);
% %             smoothed_field(:,:,:,subj) = smoothing_store;
% %         end
% %         teval_lat = mvtstat(smoothed_field, Ldim);
% %     end


% options = optimoptions(@fmincon,'Display','off'); %Ensures that no output is displayed.
% for peakI = 1:npeaks
%     peaklocs(:, peakI) = fmincon(@(tval) -tcf(tval), peak_est_locs(:, peakI), A, b, [], [], [], [], [], options);
%     peakvals(peakI) = tcf(peaklocs(:, peakI));
% end
% [ peaklocs, peakvals ] = findlms( masked_field, peak_est_locs, 1.5 );
% npeaks = size(peak_est_locs, 2);
% peaklocs = zeros(D, npeaks);
% peakvals = zeros(1, npeaks); 
% 
% A = [eye(D);-eye(D)];
% b = [Ldim(:)+0.5;ones(D,1)-0.5];
% b = zeros(2*D,1);
% for d = 1:D
%     b(d) = xvals_vecs{d}(end);
% end
% for d = 1:D
%     b(d+D) = -xvals_vecs{d}(1);
% end
% for peakI = 1:npeaks
%     peak_locs(:, peakI) = fmincon(@(tval) -tcf(tval), peak_est_locs(:, peakI), A, b);
%     peak_vals(peakI) = tcf(peak_locs(:, peakI));
% end
