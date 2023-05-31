function enlarged_xvals = expand_xvals( xvals, lower, upper )
% expand_xvals( xvals, lower, upper )
%--------------------------------------------------------------------------
% ARGUMENTS
% xvals - A cell array of vectors to be expanded.
% lower - A vector giving the number of values to add to the lower end of 
%        each vector.
% upper - A vector giving the number of values to add to the upper end of 
%         each vector.
%--------------------------------------------------------------------------
% OUTPUT
% enlarged_xvals - A cell array of the expanded vectors.
%--------------------------------------------------------------------------
% EXAMPLES
% expand_xvals({3:10, 4:5}, [1,2], [3,4])
%--------------------------------------------------------------------------
% AUTHOR: Samuel Davenport
%--------------------------------------------------------------------------
if ~exist('lower', 'var')
    lower = 1;
end
if ~exist('upper', 'var')
    upper = 1;
end

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------
enlarged_xvals = cell(1, length(xvals));
for I = 1:length(xvals)
    enlarged_xvals{I} = (xvals{I}(1)-lower(I)):(xvals{I}(end)+upper(I));
end

end

