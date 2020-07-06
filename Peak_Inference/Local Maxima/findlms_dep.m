 function [ peaklocs, peakvals ] = findlms( fn, initial_estimates, box_size, mfield )
% FINDLMS( fn, initial_estimates, box_size ) finds the local maxima of 
% function by searching within boxes centred at the initial estimates
%--------------------------------------------------------------------------
% ARGUMENTS
% fn            a function handle of a function to maximize.
% initial_estimates     a D by nestimates matrix where each column is a
%                       D-dimensional initial estimate of a local maxima location
% box_size      a 1 by nestimates cell array specifying the size of the box 
%               around each point within which to locally search.  
%               Alternately if the box_size is numeric then this box size 
%--------------------------------------------------------------------------
% OUTPUT
% peaklocs      a vector with the locations of the local maxima
%--------------------------------------------------------------------------
% EXAMPLES
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
% Calculate the number of points at which to initialize
npeaks = size(initial_estimates,2);

% Find the dimension of the data
D = size(initial_estimates,1);

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

%%  add/check optional values
%--------------------------------------------------------------------------
if ~exist('box_size', 'var')
    box_size = 1;
end

if ~iscell(box_size)
    box_size = repmat({box_size}, 1, npeaks);
end

% Need to include this bit to enable the function to find lms in non-square
% areas!
% if isequal(size(box_size), [1,1])
%     box_size = repmat(box_size, D, 1);
% elseif size(box_size,1) == 1
%     box_size = box_size'; %To ensure that it is a column vector
%     if length(box_size)~= D || size(box_size,2) ~= 1
%         error('box_size must be a column vector of length D')
%     end
% elseif length(box_size)~= D || size(box_size,2) ~= 1
%     error('box_size must be a column vector of length D')
% end


%%  main function
%--------------------------------------------------------------------------
% Set the matrix for the linear expansion.
A = [eye(D);-eye(D)];

% Initialize peak location and value matrices
peaklocs = zeros(D, npeaks);
peakvals = zeros(1, npeaks);

% Set the options for the optimization function fmincon
options = optimoptions(@fmincon,'Display','off'); % Ensures that no output 
                                                  % is displayed.
options.Algorithm = 'sqp';
% A bizarre constant that seems to be needed in order for the algorithm to 
% always converge if it is initialized at an integer!
extra = 0.00001*(floor(initial_estimates) == initial_estimates);
initial_estimates = initial_estimates + extra;

% Calculate the maximum locations on a box around each initial estimate
for peakI = 1:npeaks
    b = zeros(2*D,1);
    b(1:D) = initial_estimates(:,peakI) + box_size{peakI}' - extra(:,peakI); % Upper box limit
    b((D+1):(2*D)) = -(initial_estimates(:,peakI) - box_size{peakI}' - extra(:,peakI)); % Lower box limit
    
    % Use the optimization routine fmincon to find the peaks
    peaklocs(:,peakI) = fmincon(@(tval) -fn(tval), initial_estimates(:,peakI), A, b, [], [], [], [], mfield, options);
    
    % Evaluate the fn at the peak locations
    peakvals(peakI) = fn(peaklocs(:, peakI));
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

