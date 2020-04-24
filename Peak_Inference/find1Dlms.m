function peaklocs = find1Dlms( fn, lower_bounds, upper_bounds, initial_estimates )
% FIND1DLMS( fn, bounds ) finds the local maxima of a 1D function within
% the specified bounds starting from given initial estimates.
%--------------------------------------------------------------------------
% ARGUMENTS
% fn    a function handle of a 1D function to maximize.
% bounds    a length 2 vector: [lower_bound, upper_bound] which gives the
%           lower bound and upper bound for the search to take place.
%--------------------------------------------------------------------------
% OUTPUT
% peaklocs  a vector with the locations of the local maxima
%--------------------------------------------------------------------------
% EXAMPLES
% FWHM = 3;
% Y = [1,2,1];
% findconvpeaks(Y, FWHM, 1)
% cfield = @(tval) applyconvfield(tval, Y, 3)
% find1Dlms( cfield, [1,3], 2.5)
%
% %Boundary peaks example
% Y = [1,2,1,2,1];
% findconvpeaks(Y, FWHM, 1)
% cfield = @(tval) applyconvfield(tval, Y, 3)
% find1Dlms( cfield, 1, 5, [2.5,4.5])
% find1Dlms( cfield, [1,3], [3,5], [2.5,4.5]) %fmincon failing on the
% boundary could apply NR here to make it even more accurate, but perhaps
% we don't care.
%
% %Multiple peaks - same height
% Y = [1,2,1,1,1,1,1,2,1];
% findconvpeaks(Y, FWHM, 1)
% cfield = @(tval) applyconvfield(tval, Y, 3)
% find1Dlms( cfield, [1,length(Y)], [2.5,6.5])
% find1Dlms( cfield, [1,6], [3,9], [2.5,6.5])
%
% %Multiple peaks - different heights - still finds the local peaks!
% Y = [1,2,1,1,1,1,1,4,1];
% findconvpeaks(Y, FWHM, 1)
% cfield = @(tval) applyconvfield(tval, Y, 3)
% find1Dlms( cfield, [1,length(Y)], [2.5,6.5])
% find1Dlms( cfield, [1,6], [3,9], [2.5,6.5])
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
npeaks = length(initial_estimates);
if isequal(size(lower_bounds),[1,1])
    lower_bounds = repmat(lower_bounds, 1, npeaks);
end
if isequal(size(upper_bounds),[1,1])
    upper_bounds = repmat(upper_bounds, 1, npeaks);
end

D = 1;
A = [eye(D);-eye(D)];
peaklocs = zeros(1, npeaks);

options = optimoptions(@fmincon,'Display','off'); %Ensures that no output is displayed.
for peakI = 1:npeaks
    b = zeros(2*D,1);
    for d = 1:D
        b(d) = upper_bounds(peakI);
    end
    for d = 1:D
        b(d+D) = -lower_bounds(peakI);
    end
    peaklocs(peakI) = fmincon(@(tval) -fn(tval), initial_estimates(peakI), A, b, [], [], [], [], [], options);
end

end

