 function [maxnmin, LKCs, alphathresholds] = record_coverage( spfn, sample_size, params, niters, npeaks, version )
% RECORD_COVERAGE( data, FWHM, mask, B, sample_size ) estimates the coverage
% provided by a variety of RFT implementations including non-stationary and
% stationary convolution and lattice versions.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  spfn         a function handle which given an integer number of subjects
%               nsubj generates nsubj random fields and saves them in a
%               data type of class field
%  sample_size  the size of each sample to be sampled from the data
%  params       an object of class ConvFieldParams
% Optional
%  niters        the number of resamples of the data to do
%  npeaks       the number of peaks of the lattice data around which to
%               search for the local maxima
%  version      the version of the LKC estimation to use
%  subsets
%--------------------------------------------------------------------------
% OUTPUT
%  coverage    a structural array with entries:
%     .conv    the average coverage obtained over the number of niters
%               using the convolution approach
%     .lat     coverage obtained via evaluation on a lattice (this is
%              conservative)
%     .finelat   the coverage obtained with evaluation on the fine lattice
%              given by resadd, this is conservative but not as bad as 
%             coverage.lat. On small domains or for very smooth fields this
%             will be similar to coverage.conv, especially for high resadd.
%--------------------------------------------------------------------------
% EXAMPLES
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

% Obtain a sample of the field
single_sample_field = spfn(1);

% Obtain properties of the data
D = single_sample_field.D;

% mask = single_sample_field.mask;
% Calculate the Euler characteristic
% L0 = EulerChar(mask, 0.5, D);

%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'niters', 'var' )
    % default option of niters
    niters = 1000;
end

% Set the default do_spm value
if ~exist('version', 'var')
    if D < 3
        version = true(ones(1,D));
    elseif D == 3
        version = [true, true, false];
    end
end

if ~exist('npeaks', 'var')
    npeaks = 3;
end

% Determine whether to use fixed random subsets in advance instead of a
% given number of subjects
if iscell(sample_size)
    %Obtain the subsets
    subsets = sample_size;
    
    % Determine the sample size to use
    sample_size = length(subsets{1});
    
    % Calculate the number of subsets stored
    sample_niters = length(subsets);
    
    % Ensure that the number of iterations is the same as the number of
    % samples provided
    if sample_niters ~= niters
        error('The number of subsets must be the same as the number of iterations\n')
    end
    
    % Define an indicator that shows you're using random subsets
    use_subsets = 1;
else
    % Define an indicator that shows you're not using random subsets
    use_subsets = 0;
end

%%  Main Function Loop
%--------------------------------------------------------------------------
% Initialize vectors to store the maxima
latmaxima     = zeros( 1, niters );
finelatmaxima = zeros( 1, niters );
convmaxima    = zeros( 1, niters );

% Initialize vectors to store the minima
latminima     = zeros( 1, niters );
finelatminima = zeros( 1, niters );
convminima    = zeros( 1, niters );

% Initialize matrices to store the alpha thresholds and the LKCs
alphathresholds = zeros( 1, niters );
LKCs.L = zeros( D, niters );
LKCs.L0 = zeros(1, niters);

% Initialize matrices to store all of the high maxima and low minima
allmaxima = zeros( npeaks, niters );
allminima = zeros( npeaks, niters );

% Initialize the coverage structure
maxnmin = struct();

% Main
for b = 1:niters
    %Display b if mod(b,10) = 0
    modul(b,10)
    
    % Obtain the data
    if use_subsets
        lat_data = spfn(subsets{b});
    else
        lat_data = spfn(sample_size);
    end
    
    lat_data = Mask(lat_data);
    [ ~, threshold, maximum, L, minimum ] = vRFT(lat_data, params, npeaks, 1, version);
    LKCs.L(:,b) = L.L';
    LKCs.L0(b) = L.L0;
    if any(isnan(L))
        warning('NAN LKC recorded')
    end
    
    % Record the maxima 
    latmaxima(b) = maximum.lat;
    finelatmaxima(b) = maximum.finelat;
    convmaxima(b) = maximum.conv;
    allmaxima(1:npeaks,b) = maximum.allmaxima';
    alphathresholds(b) = threshold;
    
    % Error checking loop 
    if maximum.finelat > maximum.conv + 10^(-2)
        a = 1
    end
    
    % Record the minima
    latminima(b) = minimum.lat;
    finelatminima(b) = minimum.finelat;
    convminima(b) = minimum.conv;
    allminima(1:npeaks,b) = minimum.allmminima';
end

maxnmin.nsubj = sample_size;

% Assign the maxima
maxnmin.finelatmaxima = finelatmaxima;
maxnmin.latmaxima  = latmaxima;
maxnmin.convmaxima = convmaxima;
maxnmin.allmaxima  = allmaxima;

% Assign the minima
maxnmin.finelatminima = finelatminima;
maxnmin.latminima  = latminima;
maxnmin.convminima = convminima;
maxnmin.allminima  = allminima;

end