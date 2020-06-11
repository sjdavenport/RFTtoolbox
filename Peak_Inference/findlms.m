 function [ peaklocs, peakvals ] = findlms( fn, initial_estimates, box_size )
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
% %% 1D examples
% FWHM = 3;
% Y = [1,2,1];
% findconvpeaks(Y, FWHM, 1)
% cfield = @(tval) applyconvfield(tval, Y, 3)
% findlms( cfield, 2.5, 1 )
%
% %Multiple peaks - same height
% Y = [1,2,1,1,1,1,1,2,1]; FWMM = 3;
% findconvpeaks(Y, FWHM, 2)
% cfield = @(tval) applyconvfield(tval, Y, FWHM)
% xvals_fine = 1:0.1:length(Y);
% plot(xvals_fine, convfield(Y, FWHM, 0.1, 1))
% findlms( cfield, [2.5,6.5])
%
% %% 2D examples
% FWHM = 3; Y = [1,1,1,1;1,2,2,1;1,2,2,1;1,1,1,1];
% cfield = @(tval) applyconvfield(tval, Y, 3)
% surf(convfield(Y, FWHM, 0.1, 2))
% fine_eval = convfield(Y, 2, 0.01, 2);
% findlms( cfield, [2,2]', 1)
%
% %2D multiple peaks
% Y = [5,1,1,1;1,1,1,1;1,1,1,1;1,1,1,5]
% surf(convfield(Y, 2, 0.1, 2))
% findconvpeaks(Y, 2, [1,1;4,4]')
% cfield = @(tval) applyconvfield(tval, Y, 2)
% findlms( cfield, [1,1;4,4]', 4 )
%
% % Works with functions that take NaN values so long as the initial
% % estimate is well defined!
% mask = [0,1,1];
% mask2 = [0,1,0];
% Y = 3:-1:1;
% FWHM = 4;
% mfield = @(x) zero2nan(mask_field( x, mask2 ));
% cfield = @(x) applyconvfield(x, Y, FWHM, -1, 1:3, mask)
% masked_field = @(x) -mfield(x).*cfield(x);
% xvals_fine = 1:0.1:3;
% plot(xvals_fine, masked_field(xvals_fine))
% xlim([1,3])
% findlms( masked_field, 2, 2 )
%
% mask = [0,1,1,0,1,1];
% mask2 = [0,1,0,0,1,0];
% Y = [3:-1:1, 3:-1:1];
% nvox = length(Y);
% FWHM = 4;
% mfield = @(x) zero2nan(mask_field( x, mask2 ));
% cfield = @(x) applyconvfield(x, Y, FWHM, -1, 1:nvox, mask)
% masked_field = @(x) -mfield(x).*cfield(x);
% xvals_fine = 1:0.1:nvox;
% plot(xvals_fine, masked_field(xvals_fine))
% xlim([1,nvox])
% findlms( masked_field, 2, 2 )
% findlms( masked_field, 5, 10 ) % Note that it fails to search beyond the NaNs!
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if nargin < 3
    box_size = 1;
end

npeaks = size(initial_estimates,2);

D = size(initial_estimates,1);
try
    fn(ones(D,1));
catch
    error('The dimensions of the initial esimates do not match the dimension of the function')
end

if any(isnan(fn(initial_estimates)))
   error('All initial estimates must be well defined and lie within the mask')
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

A = [eye(D);-eye(D)];
peaklocs = zeros(D, npeaks);
peakvals = zeros(1, npeaks);

% options = optimoptions(@fmincon,'Display','off', 'Algorithm', 'sqp'); %Ensures that no output is displayed.
options = optimoptions(@fmincon,'Display','off'); %Ensures that no output is displayed.
extra = 0.0001; %Weirdly needed in order for the algorithm to always converge.
initial_estimates = initial_estimates + extra;
for peakI = 1:npeaks
    b = zeros(2*D,1);
    b(1:D) = initial_estimates(:,peakI) + box_size{peakI}';
    b((D+1):(2*D)) = -(initial_estimates(:,peakI) - box_size{peakI}');
    peaklocs(:,peakI) = fmincon(@(tval) -fn(tval), initial_estimates(:,peakI), A, b, [], [], [], [], [], options);
    peakvals(peakI) = fn(peaklocs(:, peakI));
end

end

