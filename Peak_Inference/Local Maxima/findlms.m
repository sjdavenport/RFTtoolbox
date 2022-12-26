function [ peaklocs, peakvals ] = findlms( fn, initial_estimates, ... 
                                      lowerbounds, upperbounds, algorithm )
% FINDLMS( fn, initial_estimates, box_size ) finds the local maxima of
% function by searching within boxes centred at the initial estimates
%--------------------------------------------------------------------------
% ARGUMENTS
% fn            a function handle of a function to maximize.
% initial_estimates   a D by nestimates matrix where each column is a
%                     D-dimensional initial estimate of a local maxima location
% box_size      a 1 by nestimates cell array specifying the size of the box
%               around each point within which to locally search.
%               Alternately if the box_size is numeric then this box size
%--------------------------------------------------------------------------
% OUTPUT
% peaklocs      a matrix whose columns are the locations of the local maxima
%--------------------------------------------------------------------------
% EXAMPLES
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------

if ~exist('algorithm', 'var')
   algorithm = 'sqp';
end

% Find the dimension of the data
D = size(initial_estimates,1);

% Increase the size of lowerbounds and convert to cell if necessary
if ~iscell(lowerbounds)
    if size(lowerbounds,2) < D
        lowerbounds = repmat(lowerbounds,1,D);
    end
    lowerbounds = {lowerbounds};
end
if ~iscell(upperbounds)
    if size(upperbounds,2) < D
        upperbounds = repmat(upperbounds,1,D);
    end

    upperbounds = {upperbounds};
end

% Calculate the number of points at which to initialize
npeaks = size(initial_estimates,2);

%% Error checking
%--------------------------------------------------------------------------
try
    fn(ones(D,1));
catch
    error('The dimensions of the initial estimates do not match the dimension of the function')
end

if any(isnan(fn(initial_estimates)))
    error('All initial estimates must be well defined and lie within the mask')
end

%%  main function
%--------------------------------------------------------------------------
% Initialize peak location and value matrices
peaklocs = zeros(D, npeaks);
peakvals = zeros(1, npeaks);

% Set the options for the optimization function fmincon
options = optimoptions(@fmincon,'Display','off'); % Ensures that no output

% is displayed.
options.Algorithm = algorithm;

% A bizarre constant that seems to be needed in order for the algorithm to
% always converge if it is initialized at an integer!
extra = 0.00001*(floor(initial_estimates) == initial_estimates);
initial_estimates = initial_estimates + extra;

% Calculate the maximum locations on a box around each initial estimate
for peakI = 1:npeaks
    lbs = lowerbounds{peakI};
    ubs = upperbounds{peakI};
    templocs = cell(1,size(lbs,2));
    tempvals = zeros(1,size(ubs,2));
    for J = 1:size(lbs,2) %Unsure why this second loop is necessary!
        lb = lbs(:,J);
        ub = ubs(:,J);
        
        % Use the optimization routine fmincon to find the peaks
        templocs{J} = fmincon(@(tval) -fn(tval), initial_estimates(:,peakI), [], [], [], [], lb, ub, [], options);
        
        % Evaluate the fn at the peak locations
        tempvals(J) = fn(templocs{J});
    end
    [~,maxlocind] = max(tempvals);
    peaklocs(:,peakI) = templocs{maxlocind};
    peakvals(peakI) = tempvals(maxlocind);
end

end

%  initial_estimates = initial_estimates + extras;
%
% % Calculate the maximum locations on a box around each initial estimate
% for peakI = 1:npeaks
%     b = zeros(2*D,1);
%     b(1:D) = initial_estimates(:,peakI) + box_size{peakI}' - extras(:,peakI); % Upper box limit
%     b((D+1):(2*D)) = -(initial_estimates(:,peakI) - box_size{peakI}' - extras(:,peakI)); % Lower box limit
%
%     % Use the optimization routine fmincon to find the peaks
%     peaklocs(:,peakI) = fmincon(@(tval) -fn(tval), initial_estimates(:,peakI), A, b, [], [], [], [], [], options);
%
%     % Evaluate the fn at the peak locations
%     peakvals(peakI) = fn(peaklocs(:, peakI));

% DEPRECTATED
% options = optimoptions(@fmincon,'Display','off', 'Algorithm', 'sqp'); %Ensures that no output is displayed.

