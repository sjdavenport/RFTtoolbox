function spfn = getbootspfn( data, do_gaussianize, mask )
% GETBOOTSPFN( data, mask, do_gaussianize ) obtains a function that 
% generates boostrap samples of fields from the data
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  data         an array giving the data 
% Optional
%  do_gaussianize  0/1 whether or not to Gaussianize the data. Default is 1
%                  i.e. to perform Gaussianization.
%  mask         the mask for the data. If this is an array then it is taken
%               to be the mask for the data. If is it not specified then
%               the masks are obtained from the corresponding directory
%--------------------------------------------------------------------------
% OUTPUT
%  spfn         a sample function
%--------------------------------------------------------------------------
% EXAMPLES
% %% %% 1D Examples
% %% Simple 1D example
% tot_nsubj = 30; nvox = 100; data = normrnd(0,1,nvox, tot_nsubj); 
% mask = true(1,nvox)';
% spfn = getbootspfn( data, 0, mask );
% spfn(10)
% 
% %% %% 2D Examples
% %% Simple 2D example
% tot_nsubj = 50; Dim = [20,30]; data = normrnd(0,1,[Dim,tot_nsubj]); 
% mask = true(Dim);
% spfn = getbootspfn( data, 0, mask );
% spfn(20)
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Ensure the mask is logical
mask = logical(mask);

% Calculate the dimension of the mask
D = length(size(mask));

% Cover the D = 1 case
if any(size(mask) == 1) && D == 2
    D = 1;
end

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist('do_gaussianize', 'var')
    do_gaussianize = 1;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
indexD = repmat( {':'}, 1, D );

% Obtain the size of the data
sD = size(data);

% Obtain the total number of subjects
total_nsubj = sD(end);

spfn = @(subset2use) get_sample_fields_mate(subset2use, data, do_gaussianize, mask, indexD, D, total_nsubj );

% Obtain the sample function
% if isequal(size(data), size(mask))
%     error('This is not set up')
%     spfn = @(nsubj) get_sample_fields_mate_varmask(RSfolderordata, mask, nsubj, with_rep, total_nsubj, D, indexD);
% else
% end

end

% (RSfolder, mask, nsubj, with_rep, total_nsubj, D, indexD, do_gaussianize)
% Function to obtain a subset field
function lat_data = get_sample_fields_mate(subset2use, data, do_gaussianize, mask, index, D, total_nsubj )
    % Allows for the option to choose a random (not pre-specified) set of
    % subjects with a given amount.
    if length(subset2use) == 1
        % Choose a random subset without replacement
        subset2use = randsample(total_nsubj,subset2use,0);
    end
    index{D+1} = subset2use;
    lat_data = Field(mask);
    lat_data.field = data(index{:});
    if do_gaussianize
        lat_data = Gaussianize(lat_data);
    end
end

% Function to obtain a subset field with variable masks
function out = get_sample_fields_mate_varmask(data, mask_data, nsubj, with_rep, total_nsubj, D, index)
    out.subset = randsample(total_nsubj,nsubj,with_rep);
    index{D+1} = out.subset;
    mask = logical(prod(mask_data(index{:}), D+1));
    out.lat_data = Field(mask);
    out.lat_data.field = data(index{:});
end
