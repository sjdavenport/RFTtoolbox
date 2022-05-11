function [peaklocs, peakvals, subset] = findconvpeaks(lat_data, FWHM, ...
      peak_est_locs, field_type, truncation, use_bound, use_fn_eval, meanfn, meanonlat)
% FINDCONVPEAKS( lat_data, Kernel, peak_est_locs, field_type, mask,
%                                                  truncation, xvals_vecs )
% calculates the locations of peaks in a convolution field.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  lat_data  An object of class field 
%  FWHM      The amount of smoothing used to generate the convolution field
%            (need to generalize to an arbitary kernel, i.e. to accept
%            params = ConvFieldParams)
% peak_est_locs    a D by npeaks matrix where the dth column is a vector
%                  initializing the dth peak or a cell array of length D
%                  where the dth entry is a vector containing the peaks
% field_type
% meanfn
% meanonlat
%--------------------------------------------------------------------------
% OUTPUT
% peak_locs   the true locations of the top peaks in the convolution field.
%--------------------------------------------------------------------------
% EXAMPLES
% See test_findconvpeaks_z.
%--------------------------------------------------------------------------
% TODO: Change it so that the default input into findconvpeaks is of Field
% type.
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
if ~exist('field_type', 'var')
    field_type = 'Z';
end

if ~exist('use_fn_eval', 'var')
    use_fn_eval = 1;
end
if ~exist('use_bound', 'var')
    use_bound = 0;
end

% Turn the input data into a field (with mask of ones and the standard
% xvals) if the input is not a field.
if ~isa(lat_data, 'Field')
    if length(size(lat_data)) == 2 && size(lat_data,1) == 1
        lat_data = lat_data';
    end
    lat_data = Field(lat_data, true(size(lat_data)));
end

% Allow for 'z' and 't' input
field_type = upper(field_type); 

% Ensure that field_type is of the right input
if ~ischar(field_type) || (~(strcmp(field_type, 'Z') || strcmp(field_type, 'T')))
    error('The only fields supported are mean and t fields')
end

% Obtain the number of dimensions
D = lat_data.D;

%%  add/check optional values
%--------------------------------------------------------------------------
if ~exist('peak_est_locs', 'var')
    % default option of peak_est_locs
    peak_est_locs = 1; %Default is to just consider the maximum.
end

% Ensure peak_est_locs is of the right size
if size(peak_est_locs,1) > D
    error('peak_est_locs must be of the right size')
end

if ~exist('meanfn', 'var')
    use_mean = 0;
else
    use_mean = 1;
end

%% Error checking
%--------------------------------------------------------------------------
% Ensure that the field has no NaN entries
if isnan(sum(lat_data.field(:)))
    error('Can''t yet deal with nans')
end

%%  main function
%--------------------------------------------------------------------------
% If not peak estimates are provided calculate the top ones on the lattice
if (D > 2 || isnumeric(peak_est_locs))  && (isequal(size(peak_est_locs),...
                   [1,1])) && (floor(peak_est_locs(1)) == peak_est_locs(1))
    % Set the number of maxima to define
    numberofmaxima2find = peak_est_locs;
    
    % Calculate the field on the lattice
    params = ConvFieldParams(repmat(FWHM, 1, D), 1, 0);
    if strcmp(field_type, 'Z')
        lat_eval = convfield( lat_data, params );
    elseif strcmp(field_type, 'T')
        lat_eval = convfield_t( lat_data, params );
    end
    
    % Find the top local maxima of the field on the lattice
    if use_mean
        max_indices = lmindices(lat_eval.field + meanonlat, numberofmaxima2find, lat_eval.mask);
    else
        max_indices = lmindices(lat_eval.field, numberofmaxima2find, lat_eval.mask);
    end
    
    % In D = 1 you need to transpose (CHECK THIS) Note the transpose here! 
    % It's necessary for the input to other functions.
    if D == 1
        max_indices = max_indices'; 
    end
    
    % Reset this in case there are less maxima than desired
%     numberofmaxima2find = size(max_indices,2);
    
    % Store peak locations (in terms of the actual location in the xvals)
    peak_est_locs = xvaleval(max_indices, lat_eval.xvals);
%     peak_est_locs = zeros(D, numberofmaxima2find);
%     for I = 1:D
%         peak_est_locs(I, :) = lat_data.xvals{I}(max_indices(I,:));
%     end
elseif iscell(peak_est_locs)
    % Only intended to be used in the D = 1 case
    % For D = 1, if you want to calculate multiple peaks you need to enter
    % a cell array (to differentiate between this and the top number of
    % peaks).
    peak_est_locs = cell2mat(peak_est_locs);
end

% Obtain the number of peaks
npeaks = size(peak_est_locs, 2); % Calculate the number of estimates

if use_bound == 0
    % Set field that defines the mask
    mfield = @(x) mask_field( x, lat_data.mask, 0, lat_data.xvals );
    [ lowerbounds, upperbounds ] = assign_bounds( peak_est_locs, mfield );
else
    lowerbounds = cell(1, npeaks); upperbounds = cell(1, npeaks);
    for I = 1:npeaks
        lowerbounds{I} = peak_est_locs(:,I) - use_bound;
        upperbounds{I} = peak_est_locs(:,I) + use_bound;
    end
end

% Find box around peak within which to calculate
lowerbox = zeros(D,npeaks);
upperbox = zeros(D,npeaks);
for I = 1:npeaks
    lowerbox(:,I) = min(lowerbounds{I}, [], 2);
    upperbox(:,I) = max(upperbounds{I}, [], 2);
end

% set truncation distance
if ~exist('truncation', 'var')
    truncation = 4*FWHM2sigma(FWHM);
end

%%% Note atm assumes the standard xvals_vecs!!!

% Encorporate truncation
for d = 1:D
    lowerbox(d,:) = max(lowerbox(d,:)-truncation,lat_data.xvals{d}(1));
    upperbox(d,:) = min(upperbox(d,:)+truncation,lat_data.xvals{d}(end));
end
lowerbox = floor(lowerbox);
upperbox = ceil(upperbox);

% Initialize peaklocs and vals
peaklocs = zeros(D,npeaks);
peakvals = zeros(1,npeaks);

subset = cell(1,npeaks);
% Find local maxima
for I = 1:npeaks
    % Define the local xvals_vecs
    local_xvals_vecs = cell(1,D);
    subset{I} = cell(1,D);
    for d = 1:D
        local_xvals_vecs{d} = lowerbox(d,I):upperbox(d,I);
        subset{I}{d} = local_xvals_vecs{d} - lat_data.xvals{d}(1) + 1;
    end
    
    local_mask = lat_data.mask(subset{I}{:});
    
    % For t-fields make sure to include all subjects
    if strcmp(field_type, 'T')
        subset{I}{D+1} = ':';
    end
    
    % Define local lat_data subset
    local_lat_data = lat_data(subset{I}{:});
    
    if use_fn_eval == 1
        % Define local convolution field (taking truncation = 0 as have already
        % truncated to a small box)
        if strcmp(field_type, 'Z')
            kfield = @(tval) applyconvfield(tval, local_lat_data.field, FWHM, ...
                local_mask, 0, local_xvals_vecs);
        elseif strcmp(field_type, 'T')
            kfield = @(tval) applyconvfield_t( tval, local_lat_data.field, FWHM,...
                local_mask, 0, local_xvals_vecs );
            % To include the mean here need to code it into applyconvfield_t
            % (not hard!)
        end
        if use_mean
            if strcmp(field_type, 'T')
                error('The mean hasn''t been coded for t-statistic')
            end
            cfield = @(tval) kfield(tval) + meanfn(tval);
        else
            cfield = kfield;
        end
        
        [ peakloc, peakval ] = findlms( cfield, peak_est_locs(:,I), lowerbounds{I}, upperbounds{I} );
        peakvals(I) = peakval;
        peaklocs(:,I) = peakloc;
    else
        params = ConvFieldParams(repmat(FWHM, 1, D), use_fn_eval, 0);
        if strcmp(field_type, 'T')
            localcfieldeval = convfield_t(local_lat_data, params);
        else
            localcfieldeval = convfield(local_lat_data, params);
        end
        [lm, ~, peakvals(I)] = lmindices(localcfieldeval.field, 1);
        peaklocs(:,I) = xvaleval(lm, localcfieldeval.xvals);
    end
end

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


% Obtain the boundary of the mask (if this is not supplied)
% if ~exist('boundary', 'var') || isnan(boundary)
%     boundary = bndry_voxels( logical(mask), 'full' );
% end
% 
% if ~exist('truncation', 'var')
%     truncation = -1;
% end