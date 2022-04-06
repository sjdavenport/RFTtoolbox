function obj = normalize( obj, demean, destd, fiber_d )
% STD( varargin ) does take the mean of specified fiber dimensions.
%--------------------------------------------------------------------------
% ARGUMENTS
% Mandatory
%  obj         object of type Field with fibersize > 1.
%
% Optional
%  demean      if true subtract the mean. Default true
%  destd 
%               - 0  divide by 1
%               - 1  divide by the standardeviation with factor
%                    1 / (fibersize -1). Default.
%               - 2  divide by the standardeviation with factor
%                    1 / fibersize.
%  fiber_d      positive integer denoting at which dimension in the fiber,
%               the normalization takes place.
%--------------------------------------------------------------------------
% OUTPUT
% obj  an object of class Field containing the demeaned and normalized
%      fields.
%
%--------------------------------------------------------------------------
% EXAMPLES
% % create a standard kernel object
% field = wfield( [3 3], [12 12 ] )
% % sum of last fiberdimension
% f1 = std( field )
% % sum over last  last fiberdimension different way
% f2 = std( field, 2 )
% isequal(f1, f2)
% % sum over first fiber dimension 
% f1 = std( field,1 )
%--------------------------------------------------------------------------
% Author: Fabian Telschow
%--------------------------------------------------------------------------

%%  Check mandatory input and get important constants
%--------------------------------------------------------------------------



%%  Add/check optional values
%--------------------------------------------------------------------------
if ~exist( 'demean', 'var' )
    % default option of limits
    demean = 1;
end

if ~exist( 'destd', 'var' )
    % default option of ninter
    destd = 1;
end

if ~exist( 'fiber_d', 'var' )
    % default option of ninter
    fiber_d = 1;
end



%%  Main Function Loop
%--------------------------------------------------------------------------

if demean
    obj = obj - mean(obj);
end

switch destd
    case 1
        obj = obj ./ std(obj, fiber_d, 0);
    case 2
        obj = obj ./ std(obj, fiber_d, 1);
end


end 